﻿#include "stdafx.h"
#include "CalculateForce.h"
#include "Output.h"

void Multidirect_Forcing_Method(Matrix &Fx, Matrix &Fy, Matrix &u, Matrix &v, std::list<Circle> &solidList, Param par) {

	CreateMatrix(Fx_tmp, par.N1_u, par.N2_u);
	CreateMatrix(Fy_tmp, par.N1_v, par.N2_v);

	Solids_zero_force(solidList);
	Fx = Fx * 0.;
	Fy = Fy * 0.;
	int f_max = par.N_Force;
	// correct force and velocity $f_max$ times
	for (int f = 0; f <= f_max; ++f) {
		for (auto& it : solidList) {
			it.Fr = 0.;
			it.S = 0.;
		}

		CalculateForce(Fx_tmp, Fy_tmp, solidList, u, v, par);
		Fx += Fx_tmp;
		Fy += Fy_tmp;

		//Output_U(u, "u", f, par);
		//Output_V(v, "v", f, par);
	}
}

void CalculateForce(Matrix &Fx, Matrix &Fy, std::list<Circle> &iList, Matrix& u, Matrix& v, Param par) {

	for (size_t i = 0; i < par.N1_u; ++i) {
		for (size_t j = 0; j < par.N2_u; ++j) {
			Fx[i][j] = 0.0;
		}
	}
	for (size_t i = 0; i < par.N1_v; ++i) {
		for (size_t j = 0; j < par.N2_v; ++j) {
			Fy[i][j] = 0.0;
		}
	}

	std::list<Circle>::iterator solid;

#pragma omp parallel private(solid) num_threads(1)
	{
		for (solid = iList.begin(); solid != iList.end(); solid++) {
		#pragma omp single nowait
			{
				//calculating velocities us of the solid boundary
				solid->velocities();
				solid->coordinates();
			}
		}

		for (solid = iList.begin(); solid != iList.end(); solid++) {
			#pragma omp single nowait
			{
				bool AMP = true;
				if (AMP == true)
					uf_in_Nodes    (solid->Nodes, u, v, par, solid->Nn);
				else
					uf_in_Nodes_old(solid->Nodes, u, v, par, solid->Nn);

			}
		}

		for (solid = iList.begin(); solid != iList.end(); solid++) {
			#pragma omp single nowait
			{
				for (size_t k = 0; k < solid->Nn; ++k) {
					// mass force f in Lagrange nodes
					solid->Nodes[k].f = -(solid->Nodes[k].uf - solid->Nodes[k].us) / par.d_t;
					GeomVec r = solid->Nodes[k].x / length(solid->Nodes[k].x);
					solid->Fr += dot_product(r, solid->Nodes[k].f) * solid->Nodes[k].ds;   // compression force applied to the solid
					solid->S += solid->Nodes[k].ds;
				}

				solid->Fr /= solid->S;
				for (size_t k = 0; k < solid->Nn; ++k) {
					GeomVec r = (solid->Nodes[k].x) / length(solid->Nodes[k].x);
					solid->Nodes[k].f -= solid->Fr * r;
					double dV = sqrt(par.d_x*par.d_x * solid->Nodes[k].n[1] * solid->Nodes[k].n[1]
						+ par.d_y*par.d_y * solid->Nodes[k].n[2] * solid->Nodes[k].n[2]) * solid->Nodes[k].ds;
					solid->f += solid->Nodes[k].f * dV;
					solid->tau += x_product(solid->Nodes[k].x, solid->Nodes[k].f) * dV;
				}

			}
		}

		for (solid = iList.begin(); solid != iList.end(); solid++) {
			#pragma omp single nowait
			{
				CreateMatrix(Fx_temp, par.N1_u, par.N2_u);
				CreateMatrix(Fy_temp, par.N1_v, par.N2_v);

				F_to_Euler_grid_old(solid->Nodes, Fx_temp, Fy_temp, par, solid->Nn);

				int ix_max, ix_min;
				int jx_max, jx_min;
				int iy_max, iy_min;
				int jy_max, jy_min;

				GetInfluenceArea(ix_min, ix_max, jx_min, jx_max, par.N1_u - 1, par.N2_u - 1, solid->x, int(solid->r / par.d_x) + 4, par);
				GetInfluenceArea(iy_min, iy_max, jy_min, jy_max, par.N1_v - 1, par.N2_v - 1, solid->x, int(solid->r / par.d_x) + 4, par);

				// summarizing force for Euler nodes and for solids
				for (int i = ix_min; i <= ix_max; ++i) {
					for (int j = jx_min; j <= jx_max; ++j) {
						int i_real = i_real_u(i, par);
						Fx[i_real][j] += Fx_temp[i_real][j];
					}
				}
				for (int i = iy_min; i <= iy_max; ++i) {
					for (int j = jy_min; j <= jy_max; ++j) {
						int i_real = i_real_v(i, par);
						Fy[i_real][j] += Fy_temp[i_real][j];
					}
				}
			}
		}



	}
	// copy force to non-used boundary nodes
	if (par.BC == periodical) {
		for (size_t j = 0; j < par.N2_u; ++j) {
			Fx[0         ][j] = Fx[par.N1    ][j];
			Fx[par.N1 + 2][j] = Fx[2         ][j];
			Fx[1         ][j] = Fx[par.N1 + 1][j];
		}

		for (size_t j = 0; j < par.N2_v; ++j) {
			Fy[0         ][j] = Fy[par.N1    ][j];
			Fy[par.N1 + 1][j] = Fy[1         ][j];
		}
	}

	//Output(u, u, v, Fx, Fy, par.N_step, iList, par);
	//getchar();
}

void deformation_velocity(Matrix &u, Matrix &v, Matrix &Exx, Matrix &Eyy, Matrix &Exy, Param par) {

	size_t n1 = par.N1_p;
	size_t n2 = par.N2_p;

	for (size_t i = 1; i < n1 - 1; ++i) {
		for (size_t j = 1; j < n2 - 1; ++j) {

			Exx[i][j] = (u[i + 1][j] - u[i][j]) / par.d_x;
			Eyy[i][j] = (v[i][j + 1] - v[i][j]) / par.d_y;

		}
	}

	n1 = par.N1 + 1;
	n2 = par.N2 + 1;

	for (size_t i = 0; i < n1; ++i) {
		for (size_t j = 0; j < n2; ++j) {

			Exy[i][j] = ((v[i + 1][j+1] - v[i][j+1]) / par.d_x
			           + (u[i+1][j + 1] - u[i+1][j]) / par.d_y) * 0.5;

		}
	}

}

void Solids_deformation_velocity_pressure(std::list<Circle> &Solids, Matrix &Exx, Matrix &Eyy, Matrix &Exy, Matrix &p, Param par) {

	size_t np1 = par.N1_p;
	size_t np2 = par.N2_p;

	size_t nc1 = par.N1 + 1;
	size_t nc2 = par.N2 + 1;

	for (auto& solid : Solids) {

		for (size_t k = 0; k < solid.Nn; ++k) {

			int i_max, i_min;
			int j_max, j_min;


			GetInfluenceArea(i_min, i_max, j_min, j_max, np1 - 1, np2 - 1, solid.x + solid.Nodes[k].x, 4, par);
			solid.Nodes[k].Eps(1, 1) = 0.0;
			solid.Nodes[k].Eps(2, 2) = 0.0;
			solid.Nodes[k].p         = 0.0;
			for (int i = i_min; i <= i_max; ++i) {
				for (int j = j_min; j <= j_max; ++j) {
					int i_real = i;
					if (i_real <  1      ) i_real += np1 - 2;
					if (i_real >  np1 - 2) i_real -= np1 - 2;
					GeomVec xp = x_p(i_real, j, par);
					GeomVec xs = solid.x + solid.Nodes[k].x;
					solid.Nodes[k].Eps(1, 1) += Exx[i_real][j] * DeltaFunction(xp[1] - xs[1], xp[2] - xs[2], par.d_x, par.d_y);
					solid.Nodes[k].Eps(2, 2) += Eyy[i_real][j] * DeltaFunction(xp[1] - xs[1], xp[2] - xs[2], par.d_x, par.d_y);
					solid.Nodes[k].p         +=   p[i_real][j] * DeltaFunction(xp[1] - xs[1], xp[2] - xs[2], par.d_x, par.d_y);
					if (par.BC == periodical) {
						solid.Nodes[k].Eps(1, 1) += Exx[i_real][j] * DeltaFunction(xp[1] - xs[1] - par.L, xp[2] - xs[2], par.d_x, par.d_y);
						solid.Nodes[k].Eps(1, 1) += Exx[i_real][j] * DeltaFunction(xp[1] - xs[1] + par.L, xp[2] - xs[2], par.d_x, par.d_y);
						solid.Nodes[k].Eps(2, 2) += Eyy[i_real][j] * DeltaFunction(xp[1] - xs[1] - par.L, xp[2] - xs[2], par.d_x, par.d_y);
						solid.Nodes[k].Eps(2, 2) += Eyy[i_real][j] * DeltaFunction(xp[1] - xs[1] + par.L, xp[2] - xs[2], par.d_x, par.d_y);
						solid.Nodes[k].p         +=   (p[i_real][j] + dpdx_Poiseuille(par.H, par.Re)*par.L) * DeltaFunction(xp[1] - xs[1] - par.L, xp[2] - xs[2], par.d_x, par.d_y);
						solid.Nodes[k].p         +=   (p[i_real][j] - dpdx_Poiseuille(par.H, par.Re)*par.L) * DeltaFunction(xp[1] - xs[1] + par.L, xp[2] - xs[2], par.d_x, par.d_y);
					}
				}
			}

			GetInfluenceArea(i_min, i_max, j_min, j_max, nc1 - 1, nc2 - 1, solid.x + solid.Nodes[k].x, 4, par);
			solid.Nodes[k].Eps(1, 2) = 0.0;
			for (int i = i_min; i <= i_max; ++i) {
				for (int j = j_min; j <= j_max; ++j) {
					int i_real = i;
					if (i_real <  0      ) i_real += nc1 - 1;
					if (i_real >  nc1 - 1) i_real -= nc1 - 1;
					if (i_real == 0      ) i_real  = nc1 - 1;
					GeomVec xc = x_c(i_real, j, par);
					GeomVec xs = solid.x + solid.Nodes[k].x;
					solid.Nodes[k].Eps(1, 2) += Exy[i_real][j] * DeltaFunction(xc[1] - xs[1], xc[2] - xs[2], par.d_x, par.d_y);
					if (par.BC == periodical) {
						solid.Nodes[k].Eps(1, 2) += Exy[i_real][j] * DeltaFunction(xc[1] - xs[1] - par.L, xc[2] - xs[2], par.d_x, par.d_y);
						solid.Nodes[k].Eps(1, 2) += Exy[i_real][j] * DeltaFunction(xc[1] - xs[1] + par.L, xc[2] - xs[2], par.d_x, par.d_y);
					}
				}
			}

		}
	}

}

void MatrixMultiplyAMP(std::vector<float>& vC,
	const std::vector<float>& vA,
	const std::vector<float>& vB, int M, int N, int W)
{
	array_view<const float, 2> a(M, W, vA), b(W, N, vB);
	array_view<float, 2> c(M, N, vC);
	c.discard_data();
	parallel_for_each(c.extent,
		[=](index<2> idx) restrict(amp) {
		int row = idx[0]; int col = idx[1];
		float sum = 0.0f;
		for (int i = 0; i < W; i++)
			sum += a(row, i) * b(i, col);
		c[idx] = sum;
	}
	);
}

class Node_simple
{
public:
	double x_s[4];
	double uf[4];
};



std::vector<Node_simple> Copy_Node_Simple_vector(std::vector<Node>& Nodes, int N) {

	std::vector<Node_simple> Nodes_simple(N);
	for (int k = 0; k < N; ++k) {
		Nodes_simple[k].x_s[1] = Nodes[k].x_s[1];
		Nodes_simple[k].x_s[2] = Nodes[k].x_s[2];
		Nodes_simple[k].uf[1] = Nodes[k].uf[1];
		Nodes_simple[k].uf[2] = Nodes[k].uf[2];
	}
	return Nodes_simple;
}

void Copy_Node_vector(std::vector<Node_simple>& Nodes_simple, std::vector<Node>& Nodes, int N) {

	for (int k = 0; k < N; ++k) {
		Nodes[k].x_s[1] = Nodes_simple[k].x_s[1];
		Nodes[k].x_s[2] = Nodes_simple[k].x_s[2];
		Nodes[k].uf[1] = Nodes_simple[k].uf[1];
		Nodes[k].uf[2] = Nodes_simple[k].uf[2];
	}
}

void uf_in_Nodes(std::vector<Node>& Nodes, Matrix &u, Matrix &v, Param par, int Nn)
{
	int N1_period = par.N1;
	int N1_u = par.N1_u;
	int N2_u = par.N2_u;
	int N1_v = par.N1_v;
	int N2_v = par.N2_v;
	double d_x = par.d_x;
	double d_y = par.d_y;
	boundary_conditions BC = par.BC;
	double L = par.L;
	std::clock_t begin = std::clock();

	std::vector<Node_simple> Nodes_simple = Copy_Node_Simple_vector(Nodes, Nn);
	array_view<Node_simple, 1> Nodes_AV(Nn, Nodes_simple);

	double *u_ = new double[N1_u*N2_u];
	double *v_ = new double[N1_v*N2_v];

	Matrix_to_DoubleArray(u, u_, par.BC);
	Matrix_to_DoubleArray(v, v_, par.BC);

	array_view <double, 2> u_AV(N2_u, N1_u, u_);
	array_view <double, 2> v_AV(N2_v, N1_v, v_);

	//std::cout << "Matrix u = " << u[10][20] << ";  double* = " << u_[10 + N1 * 20] << "; u_AV = " << u_AV(20, 10) << std::endl;
	//std::getchar();

	//Nodes_AV.discard_data();
	parallel_for_each(Nodes_AV.extent, [=](index<1> k) restrict(amp) {
		InfluenceArea IA_u = GetInfluenceArea_(N1_u - 1, N2_u - 1, Nodes_AV[k].x_s, 4, BC, d_x, d_y);
		//calculating fluid velocity uf in Lagrange nodes by using near Euler nodes and discrete delta function
		Nodes_AV[k].uf[1] = 0.0;
		for (int i = IA_u.i_min; i <= IA_u.i_max; ++i) {
			for (int j = IA_u.j_min; j <= IA_u.j_max; ++j) {
				int i_real = i_real_u_(i, N1_period);
				double* xu = x_u_(i, j, d_x, d_y);
				Nodes_AV[k].uf[1] += u_AV(j, i_real) * DeltaFunction_(xu[1] - Nodes_AV[k].x_s[1], xu[2] - Nodes_AV[k].x_s[2], d_x, d_y);
				if (BC == periodical) {
					Nodes_AV[k].uf[1] += u_AV(j, i_real) * DeltaFunction_(xu[1] - Nodes_AV[k].x_s[1] - L, xu[2] - Nodes_AV[k].x_s[2], d_x, d_y);
					Nodes_AV[k].uf[1] += u_AV(j, i_real) * DeltaFunction_(xu[1] - Nodes_AV[k].x_s[1] + L, xu[2] - Nodes_AV[k].x_s[2], d_x, d_y);
				}
			}
		}

		InfluenceArea IA_v = GetInfluenceArea_(N1_v - 1, N2_v - 1, Nodes_AV[k].x_s, 4, BC, d_x, d_y);
		//calculating fluid velocity uf in Lagrange nodes by using near Euler nodes and discrete delta function
		Nodes_AV[k].uf[2] = 0.0;
		for (int i = IA_v.i_min; i <= IA_v.i_max; ++i) {
			for (int j = IA_v.j_min; j <= IA_v.j_max; ++j) {
				int i_real = i_real_v_(i, N1_period);
				double* xv = x_v_(i, j, d_x, d_y);
				Nodes_AV[k].uf[2] += v_AV(j, i_real) * DeltaFunction_(xv[1] - Nodes_AV[k].x_s[1], xv[2] - Nodes_AV[k].x_s[2], d_x, d_y);
				if (BC == periodical) {
					Nodes_AV[k].uf[1] += v_AV(j, i_real) * DeltaFunction_(xv[1] - Nodes_AV[k].x_s[1] - L, xv[2] - Nodes_AV[k].x_s[2], d_x, d_y);
					Nodes_AV[k].uf[1] += v_AV(j, i_real) * DeltaFunction_(xv[1] - Nodes_AV[k].x_s[1] + L, xv[2] - Nodes_AV[k].x_s[2], d_x, d_y);
				}
			}
		}
	}
	);
	Nodes_AV.synchronize();

	std::clock_t end = std::clock();
	std::cout << "time u new " << end - begin << std::endl;

	Copy_Node_vector(Nodes_simple, Nodes, Nn);

}

void uf_in_Nodes_old(std::vector<Node>& Nodes, Matrix &u, Matrix &v, Param par, int Nn){
	
	std::clock_t begin = std::clock();

	for (size_t k = 0; k < Nn; ++k) {
		//calculating fluid velocity uf in Lagrange nodes by using near Euler nodes and discrete delta function
		int i_max, i_min;
		int j_max, j_min;
		GetInfluenceArea(i_min, i_max, j_min, j_max, par.N1_u - 1, par.N2_u - 1, Nodes[k].x_s, 4, par);

		Nodes[k].uf[1] = 0.0;
		for (int i = i_min; i <= i_max; ++i) {
			for (int j = j_min; j <= j_max; ++j) {
				int i_real = i_real_u(i, par);
				GeomVec xu = x_u(i_real, j, par);
				Nodes[k].uf[1] += u[i_real][j] * DeltaFunction(xu[1] - Nodes[k].x_s[1], xu[2] - Nodes[k].x_s[2], par.d_x, par.d_y);
				if (par.BC == periodical) {
					Nodes[k].uf[1] += u[i_real][j] * DeltaFunction(xu[1] - Nodes[k].x_s[1] - par.L, xu[2] - Nodes[k].x_s[2], par.d_x, par.d_y);
					Nodes[k].uf[1] += u[i_real][j] * DeltaFunction(xu[1] - Nodes[k].x_s[1] + par.L, xu[2] - Nodes[k].x_s[2], par.d_x, par.d_y);
				}
			}
		}
	}
	for (size_t k = 0; k < Nn; ++k) {
		//calculating fluid velocity uf in Lagrange nodes by using near Euler nodes and discrete delta function
		int i_max, i_min;
		int j_max, j_min;
		GetInfluenceArea(i_min, i_max, j_min, j_max, par.N1_v - 1, par.N2_v - 1, Nodes[k].x_s, 4, par);
		Nodes[k].uf[2] = 0.0;
		for (int i = i_min; i <= i_max; ++i) {
			for (int j = j_min; j <= j_max; ++j) {
				int i_real = i_real_v(i, par);
				GeomVec xv = x_v(i_real, j, par);
				Nodes[k].uf[2] += v[i_real][j] * DeltaFunction(xv[1] - Nodes[k].x_s[1], xv[2] - Nodes[k].x_s[2], par.d_x, par.d_y);
				if (par.BC == periodical) {
					Nodes[k].uf[2] += v[i_real][j] * DeltaFunction(xv[1] - Nodes[k].x_s[1] - par.L, xv[2] - Nodes[k].x_s[2], par.d_x, par.d_y);
					Nodes[k].uf[2] += v[i_real][j] * DeltaFunction(xv[1] - Nodes[k].x_s[1] + par.L, xv[2] - Nodes[k].x_s[2], par.d_x, par.d_y);
				}
			}
		}
	}
	std::clock_t end = std::clock();
	std::cout << "time u old " << end - begin << std::endl;
}


void F_to_Euler_grid(std::vector<Node>& Nodes, Matrix &Fx_temp, Matrix &Fy_temp, Param par, int Nn) {

	int N1_u = par.N1_u;
	int N2_u = par.N2_u;
	int N1_v = par.N1_v;
	int N2_v = par.N2_v;

	std::clock_t begin = std::clock();

	std::vector<Node_simple> Nodes_simple = Copy_Node_Simple_vector(Nodes, Nn);
	array_view<Node_simple, 1> Nodes_AV(Nn, Nodes_simple);

	double *Fx_temp_ = new double[N1_u*N2_u];
	double *Fy_temp_ = new double[N1_v*N2_v];

	Matrix_to_DoubleArray(Fx_temp, Fx_temp_, par.BC);
	Matrix_to_DoubleArray(Fy_temp, Fy_temp_, par.BC);

	array_view <double, 2> u_AV(N2_u, N1_u, Fx_temp_);
	array_view <double, 2> v_AV(N2_v, N1_v, Fy_temp_);

	for (size_t k = 0; k < Nn; ++k) {

		int ix_max, ix_min;
		int jx_max, jx_min;

		int iy_max, iy_min;
		int jy_max, jy_min;

		GetInfluenceArea(ix_min, ix_max, jx_min, jx_max, par.N1_u - 1, par.N2_u - 1, Nodes[k].x_s, 4, par);
		GetInfluenceArea(iy_min, iy_max, jy_min, jy_max, par.N1_v - 1, par.N2_v - 1, Nodes[k].x_s, 4, par);

		//calculating velocities us of the solid boundary
		double dn = sqrt(par.d_x*par.d_x * Nodes[k].n[1] * Nodes[k].n[1]
		               + par.d_y*par.d_y * Nodes[k].n[2] * Nodes[k].n[2]) / par.d_x / par.d_y;
		// calculating force force_temp for Euler nodes caused by k-th solid
		for (int i = ix_min; i <= ix_max; ++i) {
			for (int j = jx_min; j <= jx_max; ++j) {
				int i_real = i_real_u(i, par);
				GeomVec xu = x_u(i_real, j, par);
				Fx_temp[i_real][j] += Nodes[k].f[1] * DeltaFunction(xu[1] - Nodes[k].x_s[1], xu[2] - Nodes[k].x_s[2], par.d_x, par.d_y) * dn * Nodes[k].ds;

				if (par.BC == periodical) {

					GeomVec xu_plus = xu;
					xu_plus[1] += par.L;
					Fx_temp[i_real][j] += Nodes[k].f[1] * DeltaFunction(xu_plus[1] - Nodes[k].x_s[1], xu_plus[2] - Nodes[k].x_s[2], par.d_x, par.d_y) * dn * Nodes[k].ds;

					GeomVec xu_minus = xu;
					xu_minus[1] -= par.L;
					Fx_temp[i_real][j] += Nodes[k].f[1] * DeltaFunction(xu_minus[1] - Nodes[k].x_s[1], xu_minus[2] - Nodes[k].x_s[2], par.d_x, par.d_y) * dn * Nodes[k].ds;

				}
			}
		}

		for (int i = iy_min; i <= iy_max; ++i) {
			for (int j = jy_min; j <= jy_max; ++j) {
				int i_real = i_real_v(i, par);
				GeomVec xv = x_v(i_real, j, par);
				Fy_temp[i_real][j] += Nodes[k].f[2] * DeltaFunction(xv[1] - Nodes[k].x_s[1], xv[2] - Nodes[k].x_s[2], par.d_x, par.d_y) * dn * Nodes[k].ds;

				if (par.BC == periodical) {

					GeomVec xv_plus = xv;
					xv_plus[1] += par.L;
					Fy_temp[i_real][j] += Nodes[k].f[2] * DeltaFunction(xv_plus[1] - Nodes[k].x_s[1], xv_plus[2] - Nodes[k].x_s[2], par.d_x, par.d_y) * dn * Nodes[k].ds;

					GeomVec xv_minus = xv;
					xv_minus[1] -= par.L;
					Fy_temp[i_real][j] += Nodes[k].f[2] * DeltaFunction(xv_minus[1] - Nodes[k].x_s[1], xv_minus[2] - Nodes[k].x_s[2], par.d_x, par.d_y) * dn * Nodes[k].ds;

				}
			}
		}
	}
	std::clock_t end = std::clock();
	std::cout << "time force " << end - begin << std::endl;
}

void F_to_Euler_grid_old(std::vector<Node>& Nodes, Matrix &Fx_temp, Matrix &Fy_temp, Param par, int Nn) {

	std::clock_t begin = std::clock();

	for (size_t k = 0; k < Nn; ++k) {

		int ix_max, ix_min;
		int jx_max, jx_min;

		int iy_max, iy_min;
		int jy_max, jy_min;

		GetInfluenceArea(ix_min, ix_max, jx_min, jx_max, par.N1_u - 1, par.N2_u - 1, Nodes[k].x_s, 4, par);
		GetInfluenceArea(iy_min, iy_max, jy_min, jy_max, par.N1_v - 1, par.N2_v - 1, Nodes[k].x_s, 4, par);

		//calculating velocities us of the solid boundary
		double dn = sqrt(par.d_x*par.d_x * Nodes[k].n[1] * Nodes[k].n[1]
		               + par.d_y*par.d_y * Nodes[k].n[2] * Nodes[k].n[2]) / par.d_x / par.d_y;
		// calculating force force_temp for Euler nodes caused by k-th solid
		for (int i = ix_min; i <= ix_max; ++i) {
			for (int j = jx_min; j <= jx_max; ++j) {
				int i_real = i_real_u(i, par);
				GeomVec xu = x_u(i_real, j, par);
				Fx_temp[i_real][j] += Nodes[k].f[1] * DeltaFunction(xu[1] - Nodes[k].x_s[1], xu[2] - Nodes[k].x_s[2], par.d_x, par.d_y) * dn * Nodes[k].ds;

				if (par.BC == periodical) {

					GeomVec xu_plus = xu;
					xu_plus[1] += par.L;
					Fx_temp[i_real][j] += Nodes[k].f[1] * DeltaFunction(xu_plus[1] - Nodes[k].x_s[1], xu_plus[2] - Nodes[k].x_s[2], par.d_x, par.d_y) * dn * Nodes[k].ds;

					GeomVec xu_minus = xu;
					xu_minus[1] -= par.L;
					Fx_temp[i_real][j] += Nodes[k].f[1] * DeltaFunction(xu_minus[1] - Nodes[k].x_s[1], xu_minus[2] - Nodes[k].x_s[2], par.d_x, par.d_y) * dn * Nodes[k].ds;

				}
			}
		}

		for (int i = iy_min; i <= iy_max; ++i) {
			for (int j = jy_min; j <= jy_max; ++j) {
				int i_real = i_real_v(i, par);
				GeomVec xv = x_v(i_real, j, par);
				Fy_temp[i_real][j] += Nodes[k].f[2] * DeltaFunction(xv[1] - Nodes[k].x_s[1], xv[2] - Nodes[k].x_s[2], par.d_x, par.d_y) * dn * Nodes[k].ds;

				if (par.BC == periodical) {

					GeomVec xv_plus = xv;
					xv_plus[1] += par.L;
					Fy_temp[i_real][j] += Nodes[k].f[2] * DeltaFunction(xv_plus[1] - Nodes[k].x_s[1], xv_plus[2] - Nodes[k].x_s[2], par.d_x, par.d_y) * dn * Nodes[k].ds;

					GeomVec xv_minus = xv;
					xv_minus[1] -= par.L;
					Fy_temp[i_real][j] += Nodes[k].f[2] * DeltaFunction(xv_minus[1] - Nodes[k].x_s[1], xv_minus[2] - Nodes[k].x_s[2], par.d_x, par.d_y) * dn * Nodes[k].ds;

				}
			}
		}
	}
	std::clock_t end = std::clock();
	std::cout << "time force " << end - begin << std::endl;
}

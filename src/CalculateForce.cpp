#include "stdafx.h"
#include "CalculateForce.h"
#include "Output.h"


void CalculateForce(Matrix &dFx, Matrix &dFy, int* Ax_beg, int* Ax_end, int* Ay_beg, int* Ay_end, std::vector<Circle> &iList, std::vector<Node> &Nodes, Matrix& u, Matrix& v, Param par) {

	int num_thr = omp_get_max_threads();

	std::clock_t begin, end;

	std::vector<Circle>::iterator solid;

	#pragma omp parallel private(solid) num_threads(num_thr)
	{
		for (solid = iList.begin(); solid != iList.end(); solid++) {
		#pragma omp single nowait
			{
				//calculating velocities us of the solid boundary
				velocities(solid,Nodes);
				coordinates(solid,Nodes);
			}
		}
	}

	//if (par.AMP == true) {
		begin = std::clock();
			uf_in_Nodes(Nodes, u, v, par, par.Nn_max);
		end = std::clock();
		std::cout << "time u new " << end - begin << std::endl;
	//}

	//if (par.AMP == false) {
		begin = std::clock();
			uf_in_Nodes_old(Nodes, u, v, par, par.Nn_max);
		end = std::clock();
		std::cout << "time u old " << end - begin << std::endl;
	//}

    #pragma omp parallel private(solid) num_threads(num_thr)
	{

		for (solid = iList.begin(); solid != iList.end(); solid++) {
		#pragma omp single nowait
			{
				for (size_t k = 0; k < solid->Nn; ++k) {
					// mass force f in Lagrange nodes
					int Ind = solid->IndNodes[k];
					Nodes[Ind].f = -(Nodes[Ind].uf - Nodes[Ind].us) / par.d_t;
					GeomVec r = Nodes[Ind].x / length(Nodes[Ind].x);
					solid->Fr += dot_product(r, Nodes[Ind].f) * Nodes[Ind].ds;   // compression force applied to the solid
					solid->S += Nodes[Ind].ds;
				}

				solid->Fr /= solid->S;
				for (size_t k = 0; k < solid->Nn; ++k) {
					int Ind = solid->IndNodes[k];
					GeomVec r = (Nodes[Ind].x) / length(Nodes[Ind].x);
					Nodes[Ind].f -= solid->Fr * r;
					solid->f += Nodes[Ind].f * Nodes[Ind].ds;
					solid->tau += x_product(Nodes[Ind].x, Nodes[Ind].f) * Nodes[Ind].ds;
				}

			}
		}
	}

	for (size_t i = 0; i < par.N1_u; ++i) {
		for (size_t j = 0; j < par.N2_u; ++j) {
			dFx[i][j] = 0.0;
		}
	}
	for (size_t i = 0; i < par.N1_v; ++i) {
		for (size_t j = 0; j < par.N2_v; ++j) {
			dFy[i][j] = 0.0;
		}
	}

	//if (par.AMP == true) {
		begin = std::clock();
			F_to_Euler_grid(Nodes, dFx, dFy, Ax_beg, Ax_end, Ay_beg, Ay_end, par, par.Nn_max);
		end = std::clock();
		std::cout << "time f new " << end - begin << std::endl;
	//}

	//if (par.AMP == false) {
		begin = std::clock();
			F_to_Euler_grid_old(Nodes, dFx, dFy, par, par.Nn_max);
		end = std::clock();
		std::cout << "time f old " << end - begin << std::endl;
	//}


	// copy force to non-used boundary nodes
	if (par.BC == periodical) {
		for (size_t j = 0; j < par.N2_u; ++j) {
			dFx[0         ][j] = dFx[par.N1    ][j];
			dFx[par.N1 + 2][j] = dFx[2         ][j];
			dFx[1         ][j] = dFx[par.N1 + 1][j];
		}

		for (size_t j = 0; j < par.N2_v; ++j) {
			dFy[0         ][j] = dFy[par.N1    ][j];
			dFy[par.N1 + 1][j] = dFy[1         ][j];
		}
	}

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

void Solids_deformation_velocity_pressure(std::vector<Circle> &Solids, std::vector<Node> &Nodes, Matrix &Exx, Matrix &Eyy, Matrix &Exy, Matrix &p, Param par) {

	size_t np1 = par.N1_p;
	size_t np2 = par.N2_p;

	size_t nc1 = par.N1 + 1;
	size_t nc2 = par.N2 + 1;

	for (auto& solid : Solids) {

		for (size_t k = 0; k < solid.Nn; ++k) {

			int i_max, i_min;
			int j_max, j_min;
			int Ind = solid.IndNodes[k];

			GetInfluenceArea(i_min, i_max, j_min, j_max, np1 - 1, np2 - 1, solid.x + Nodes[Ind].x, 4, par);
			Nodes[Ind].Eps(1, 1) = 0.0;
			Nodes[Ind].Eps(2, 2) = 0.0;
			Nodes[Ind].p         = 0.0;
			for (int i = i_min; i <= i_max; ++i) {
				for (int j = j_min; j <= j_max; ++j) {
					int i_real = i;
					if (i_real <  1      ) i_real += np1 - 2;
					if (i_real >  np1 - 2) i_real -= np1 - 2;
					GeomVec xp = x_p(i_real, j, par);
					GeomVec xs = solid.x + Nodes[Ind].x;
					Nodes[Ind].Eps(1, 1) += Exx[i_real][j] * DeltaFunction(xp[1] - xs[1], xp[2] - xs[2], par.d_x, par.d_y);
					Nodes[Ind].Eps(2, 2) += Eyy[i_real][j] * DeltaFunction(xp[1] - xs[1], xp[2] - xs[2], par.d_x, par.d_y);
					Nodes[Ind].p         +=   p[i_real][j] * DeltaFunction(xp[1] - xs[1], xp[2] - xs[2], par.d_x, par.d_y);
					if (par.BC == periodical) {
						Nodes[Ind].Eps(1, 1) += Exx[i_real][j] * DeltaFunction(xp[1] - xs[1] - par.L, xp[2] - xs[2], par.d_x, par.d_y);
						Nodes[Ind].Eps(1, 1) += Exx[i_real][j] * DeltaFunction(xp[1] - xs[1] + par.L, xp[2] - xs[2], par.d_x, par.d_y);
						Nodes[Ind].Eps(2, 2) += Eyy[i_real][j] * DeltaFunction(xp[1] - xs[1] - par.L, xp[2] - xs[2], par.d_x, par.d_y);
						Nodes[Ind].Eps(2, 2) += Eyy[i_real][j] * DeltaFunction(xp[1] - xs[1] + par.L, xp[2] - xs[2], par.d_x, par.d_y);
						Nodes[Ind].p         +=   (p[i_real][j] + dpdx_Poiseuille(par.H, par.Re)*par.L) * DeltaFunction(xp[1] - xs[1] - par.L, xp[2] - xs[2], par.d_x, par.d_y);
						Nodes[Ind].p         +=   (p[i_real][j] - dpdx_Poiseuille(par.H, par.Re)*par.L) * DeltaFunction(xp[1] - xs[1] + par.L, xp[2] - xs[2], par.d_x, par.d_y);
					}
				}
			}

			GetInfluenceArea(i_min, i_max, j_min, j_max, nc1 - 1, nc2 - 1, solid.x + Nodes[Ind].x, 4, par);
			Nodes[Ind].Eps(1, 2) = 0.0;
			for (int i = i_min; i <= i_max; ++i) {
				for (int j = j_min; j <= j_max; ++j) {
					int i_real = i;
					if (i_real <  0      ) i_real += nc1 - 1;
					if (i_real >  nc1 - 1) i_real -= nc1 - 1;
					if (i_real == 0      ) i_real  = nc1 - 1;
					GeomVec xc = x_c(i_real, j, par);
					GeomVec xs = solid.x + Nodes[Ind].x;
					Nodes[Ind].Eps(1, 2) += Exy[i_real][j] * DeltaFunction(xc[1] - xs[1], xc[2] - xs[2], par.d_x, par.d_y);
					if (par.BC == periodical) {
						Nodes[Ind].Eps(1, 2) += Exy[i_real][j] * DeltaFunction(xc[1] - xs[1] - par.L, xc[2] - xs[2], par.d_x, par.d_y);
						Nodes[Ind].Eps(1, 2) += Exy[i_real][j] * DeltaFunction(xc[1] - xs[1] + par.L, xc[2] - xs[2], par.d_x, par.d_y);
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
	double f[4];

	double ds;
};



std::vector<Node_simple> Copy_Node_Simple_vector(std::vector<Node>& Nodes, int N1, int N2) {

	std::vector<Node_simple> Nodes_simple(N2-N1);
	for (int k = 0; k < N2-N1; ++k) {
		Nodes_simple[k].x_s[1] = Nodes[N1 +k ].x_s[1];
		Nodes_simple[k].x_s[2] = Nodes[N1 + k].x_s[2];
		Nodes_simple[k].uf[1] = Nodes[N1 + k].uf[1];
		Nodes_simple[k].uf[2] = Nodes[N1 + k].uf[2];
		Nodes_simple[k].f[1] = Nodes[N1 + k].f[1];
		Nodes_simple[k].f[2] = Nodes[N1 + k].f[2];

		Nodes_simple[k].ds = Nodes[N1 + k].ds;
	}
	return Nodes_simple;
}

void Copy_Node_vector(std::vector<Node_simple>& Nodes_simple, std::vector<Node>& Nodes, int N1, int N2) {

	for (size_t k = 0; k < N2-N1; ++k) {
		Nodes[N1 + k].x_s[1] = Nodes_simple[k].x_s[1];
		Nodes[N1 + k].x_s[2] = Nodes_simple[k].x_s[2];
		Nodes[N1 + k].uf[1] = Nodes_simple[k].uf[1];
		Nodes[N1 + k].uf[2] = Nodes_simple[k].uf[2];
		Nodes[N1 + k].f[1] = Nodes_simple[k].f[1];
		Nodes[N1 + k].f[2] = Nodes_simple[k].f[2];
		Nodes[N1 + k].ds = Nodes_simple[k].ds;
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
	//std::clock_t begin = std::clock();

	std::vector<Node_simple> Nodes_simple = Copy_Node_Simple_vector(Nodes, 0, Nn);
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

	//std::clock_t end = std::clock();
	//std::cout << "time u new " << end - begin << std::endl;

	Copy_Node_vector(Nodes_simple, Nodes, 0, Nn);

}

void uf_in_Nodes_old(std::vector<Node>& Nodes, Matrix &u, Matrix &v, Param par, int Nn){
	
	//std::clock_t begin = std::clock();

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
	//std::clock_t end = std::clock();
	//std::cout << "time u old " << end - begin << std::endl;
}

void Make_interaction_Matrix(int* A_beg, int* A_end, int N1, int N2, double d_x, double d_y, std::vector<Node>& Nodes, int Nn_max, Direction Dir) {

	array_view <int, 2> A_beg_(N2, N1, A_beg);
	array_view <int, 2> A_end_(N2, N1, A_end);

	std::vector<Node_simple> Nodes_simple = Copy_Node_Simple_vector(Nodes, 0, Nn_max);
	array_view<Node_simple, 1> Nodes_AV(Nn_max, Nodes_simple);

	parallel_for_each(A_beg_.extent, [=](index<2> idx) restrict(amp) {
		int i = idx[1];
		int j = idx[0];

		double* x;
		int i_real = i;

		A_beg_(j, i_real) = 0;
		A_end_(j, i_real) = 0;

		if (Dir == Du) {
			//i_real = i_real_u_(i, N1_period);
			x = x_u_(i_real, j, d_x, d_y);
		}
		else if (Dir == Dv) {
			//i_real = i_real_v_(i, N1_period);
			x = x_v_(i_real, j, d_x, d_y);
		}

		for (int k = 0; k < Nn_max; ++k) {
			if (fabs(x[1] - Nodes_AV[k].x_s[1]) < 2.1*d_x && fabs(x[2] - Nodes_AV[k].x_s[2]) < 2.1*d_y) {
				A_beg_(j, i_real) = k;
				break;
				/*if (BC == periodical) {

					double* xu_plus = xu;
					xu_plus[1] += L;
					A_beg_(j, i_real) = k;

					double* xu_minus = xu;
					xu_minus[1] -= L;
					A_beg_(j, i_real) = k;
				}*/
			}
		}

		for (int k = Nn_max-1; k >= 0; --k) {
			if (fabs(x[1] - Nodes_AV[k].x_s[1]) < 2.1*d_x && fabs(x[2] - Nodes_AV[k].x_s[2]) < 2.1*d_y) {
				A_end_(j, i_real) = k;
				break;
				/*if (BC == periodical) {

				double* xu_plus = xu;
				xu_plus[1] += L;
				A_end_(j, i_real) = k;

				double* xu_minus = xu;
				xu_minus[1] -= L;
				A_end_(j, i_real) = k;
				}*/
			}
		}
	}
	);
}

void F_to_Euler_grid(std::vector<Node>& Nodes, Matrix &Fx_temp, Matrix &Fy_temp, int* Ax_beg, int* Ax_end, int* Ay_beg, int* Ay_end, Param par, int Nn) {

	int N1_period = par.N1;
	int N1_u = par.N1_u;
	int N2_u = par.N2_u;
	int N1_v = par.N1_v;
	int N2_v = par.N2_v;
	double d_x = par.d_x;
	double d_y = par.d_y;
	boundary_conditions BC = par.BC;
	double L = par.L;
	double l_dxdy = 1. / d_x / d_y;

	std::clock_t begin, end;

	array_view <int, 2> Ax_beg_(N2_u, N1_u, Ax_beg);
	array_view <int, 2> Ax_end_(N2_u, N1_u, Ax_end);
	array_view <int, 2> Ay_beg_(N2_v, N1_v, Ay_beg);
	array_view <int, 2> Ay_end_(N2_v, N1_v, Ay_end);

	double *Fx_temp_ = new double[N1_u*N2_u];
	double *Fy_temp_ = new double[N1_v*N2_v];

	Matrix_to_DoubleArray(Fx_temp, Fx_temp_, par.BC);
	Matrix_to_DoubleArray(Fy_temp, Fy_temp_, par.BC);

	array_view <double, 2> Fx_temp_AV(N2_u, N1_u, Fx_temp_);
	array_view <double, 2> Fy_temp_AV(N2_v, N1_v, Fy_temp_);

	std::vector<Node_simple> Nodes_simple = Copy_Node_Simple_vector(Nodes, 0, Nn);
	array_view<Node_simple, 1> Nodes_AV(Nn, Nodes_simple);

		//Fx_temp_AV.discard_data();
		parallel_for_each(Fx_temp_AV.extent, [=](index<2> idx) restrict(amp) {
			int i = idx[1];
			int j = idx[0];
			int i_real = i_real_u_(i, N1_period);
			double* xu = x_u_(i, j, d_x, d_y);

			for (int k = Ax_beg_(j, i_real); k <= Ax_end_(j, i_real); ++k) {
				if (fabs(xu[1] - Nodes_AV[k].x_s[1]) < 2.1*d_x && fabs(xu[2] - Nodes_AV[k].x_s[2]) < 2.1*d_y) {
					Fx_temp_AV(j, i_real) += Nodes_AV[k].f[1] * DeltaFunction_(xu[1] - Nodes_AV[k].x_s[1], xu[2] - Nodes_AV[k].x_s[2], d_x, d_y) * l_dxdy * Nodes_AV[k].ds;
				}
					/*if (BC == periodical) {

						double* xu_plus = xu;
						xu_plus[1] += L;
						Fx_temp_AV(j, i_real) += Nodes_AV[k].f[1] * DeltaFunction_(xu_plus[1] - Nodes_AV[k].x_s[1], xu_plus[2] - Nodes_AV[k].x_s[2], d_x, d_y) * l_dxdy * Nodes_AV[k].ds;

						double* xu_minus = xu;
						xu_minus[1] -= L;
						Fx_temp_AV(j, i_real) += Nodes_AV[k].f[1] * DeltaFunction_(xu_minus[1] - Nodes_AV[k].x_s[1], xu_minus[2] - Nodes_AV[k].x_s[2], d_x, d_y) * l_dxdy * Nodes_AV[k].ds;
					}*/
			}
		}
		);

		//Fy_temp_AV.discard_data();
		parallel_for_each(Fy_temp_AV.extent, [=](index<2> idx) restrict(amp) {
			int i = idx[1];
			int j = idx[0];
			int i_real = i_real_v_(i, N1_period);
			double* xv = x_v_(i_real, j, d_x, d_y);

			for (int k = Ay_beg_(j, i_real); k <= Ay_end_(j, i_real); ++k) {
				if (fabs(xv[1] - Nodes_AV[k].x_s[1]) < 2.1*d_x && fabs(xv[2] - Nodes_AV[k].x_s[2]) < 2.1*d_y) {
					Fy_temp_AV(j, i_real) += Nodes_AV[k].f[2] * DeltaFunction_(xv[1] - Nodes_AV[k].x_s[1], xv[2] - Nodes_AV[k].x_s[2], d_x, d_y) * l_dxdy * Nodes_AV[k].ds;
				}
					/*if (BC == periodical) {

						double* xv_plus = xv;
						xv_plus[1] += L;
						Fy_temp_AV(j, i_real) += Nodes_AV[k].f[2] * DeltaFunction_(xv_plus[1] - Nodes_AV[k].x_s[1], xv_plus[2] - Nodes_AV[k].x_s[2], d_x, d_y) * l_dxdy * Nodes_AV[k].ds;

						double* xv_minus = xv;
						xv_minus[1] -= L;
						Fy_temp_AV(j, i_real) += Nodes_AV[k].f[2] * DeltaFunction_(xv_minus[1] - Nodes_AV[k].x_s[1], xv_minus[2] - Nodes_AV[k].x_s[2], d_x, d_y) * l_dxdy * Nodes_AV[k].ds;
					}*/
			}
		}
		);

		begin = std::clock();
		Fx_temp_AV.synchronize();
		end = std::clock();
		//std::cout << "time copy Fx" << end - begin << std::endl;

		DoubleArray_to_Matrix(Fx_temp_, Fx_temp, par.BC);

		//Output_2DArray(Fx_temp_, N1_u, N2_u, "Result/", "Array", 555);
		//std::getchar();

		begin = std::clock();
		Fy_temp_AV.synchronize();
		end = std::clock();
		//std::cout << "time copy Fy" << end - begin << std::endl;

		DoubleArray_to_Matrix(Fy_temp_, Fy_temp, par.BC);

}

void F_to_Euler_grid_old(std::vector<Node>& Nodes, Matrix &Fx_temp, Matrix &Fy_temp, Param par, int Nn) {

	//std::clock_t begin = std::clock();

	//for (size_t i = 0; i < par.N1_u; ++i) {
	//	for (size_t j = 0; j < par.N2_u; ++j) {
	//		Fx_temp[i][j] = 0.0;
	//	}
	//}
	//for (size_t i = 0; i < par.N1_v; ++i) {
	//	for (size_t j = 0; j < par.N2_v; ++j) {
	//		Fy_temp[i][j] = 0.0;
	//	}
	//}

	for (size_t k = 0; k < Nn; ++k) {

		int ix_max, ix_min;
		int jx_max, jx_min;

		int iy_max, iy_min;
		int jy_max, jy_min;

		GetInfluenceArea(ix_min, ix_max, jx_min, jx_max, par.N1_u - 1, par.N2_u - 1, Nodes[k].x_s, 4, par);
		GetInfluenceArea(iy_min, iy_max, jy_min, jy_max, par.N1_v - 1, par.N2_v - 1, Nodes[k].x_s, 4, par);

		//calculating velocities us of the solid boundary
		double l_dxdy = 1. / par.d_x / par.d_y;
		// calculating force force_temp for Euler nodes caused by k-th solid
		for (int i = ix_min; i <= ix_max; ++i) {
			for (int j = jx_min; j <= jx_max; ++j) {
				int i_real = i_real_u(i, par);
				GeomVec xu = x_u(i_real, j, par);
				Fx_temp[i_real][j] += Nodes[k].f[1] * DeltaFunction(xu[1] - Nodes[k].x_s[1], xu[2] - Nodes[k].x_s[2], par.d_x, par.d_y) * l_dxdy * Nodes[k].ds;

				if (par.BC == periodical) {

					GeomVec xu_plus = xu;
					xu_plus[1] += par.L;
					Fx_temp[i_real][j] += Nodes[k].f[1] * DeltaFunction(xu_plus[1] - Nodes[k].x_s[1], xu_plus[2] - Nodes[k].x_s[2], par.d_x, par.d_y) * l_dxdy * Nodes[k].ds;

					GeomVec xu_minus = xu;
					xu_minus[1] -= par.L;
					Fx_temp[i_real][j] += Nodes[k].f[1] * DeltaFunction(xu_minus[1] - Nodes[k].x_s[1], xu_minus[2] - Nodes[k].x_s[2], par.d_x, par.d_y) * l_dxdy * Nodes[k].ds;

				}
			}
		}

		for (int i = iy_min; i <= iy_max; ++i) {
			for (int j = jy_min; j <= jy_max; ++j) {
				int i_real = i_real_v(i, par);
				GeomVec xv = x_v(i_real, j, par);
				Fy_temp[i_real][j] += Nodes[k].f[2] * DeltaFunction(xv[1] - Nodes[k].x_s[1], xv[2] - Nodes[k].x_s[2], par.d_x, par.d_y) * l_dxdy * Nodes[k].ds;

				if (par.BC == periodical) {

					GeomVec xv_plus = xv;
					xv_plus[1] += par.L;
					Fy_temp[i_real][j] += Nodes[k].f[2] * DeltaFunction(xv_plus[1] - Nodes[k].x_s[1], xv_plus[2] - Nodes[k].x_s[2], par.d_x, par.d_y) * l_dxdy * Nodes[k].ds;

					GeomVec xv_minus = xv;
					xv_minus[1] -= par.L;
					Fy_temp[i_real][j] += Nodes[k].f[2] * DeltaFunction(xv_minus[1] - Nodes[k].x_s[1], xv_minus[2] - Nodes[k].x_s[2], par.d_x, par.d_y) * l_dxdy * Nodes[k].ds;

				}
			}
		}
	}
	//std::clock_t end = std::clock();
	//std::cout << "time f old " << end - begin << std::endl;
}

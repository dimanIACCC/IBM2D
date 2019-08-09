#include "stdafx.h"
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
		u += Fx_tmp * (par.d_t);
		v += Fy_tmp * (par.d_t);
		Fx += Fx_tmp;
		Fy += Fy_tmp;

		//Output_U(u, "u", f, par);
		//Output_V(v, "v", f, par);
	}
}

void CalculateForce(Matrix &Fx, Matrix &Fy, std::list<Circle> &iList, Matrix& u, Matrix& v, Param par) {

	double dx3 = 3 * par.d_x;
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

#pragma omp parallel private(solid) num_threads(4)
	{
		for (solid = iList.begin(); solid != iList.end(); solid++) {
		#pragma omp single nowait
			{
				//calculating velocities us of the solid boundary
				solid->velocities();

				std::clock_t begin = std::clock();
				for (size_t k = 0; k < solid->Nn; ++k) {

					int i_max, i_min;
					int j_max, j_min;
					GeomVec xs = solid->x + solid->Nodes[k].x;

					//calculating fluid velocity uf in Lagrange nodes by using near Euler nodes and discrete delta function
					GetInfluenceArea(i_min, i_max, j_min, j_max, par.N1_u - 1, par.N2_u - 1, xs, 4, par);
					solid->Nodes[k].uf[1] = 0.0;
					for (int i = i_min; i <= i_max; ++i) {
						int i_real = i_real_u(i, par);
						for (int j = j_min; j <= j_max; ++j) {
							GeomVec r = x_u(i_real, j, par) - xs;
							solid->Nodes[k].uf[1] += u[i_real][j] * DeltaFunction(r[1], r[2], par);
							if (par.BC == periodical) {
								double x_minus = r[1] - par.L;
								double x_plus  = r[1] + par.L;
								if (x_minus < dx3) solid->Nodes[k].uf[1] += u[i_real][j] * DeltaFunction(x_minus, r[2], par);
								if (x_plus  < dx3) solid->Nodes[k].uf[1] += u[i_real][j] * DeltaFunction(x_plus , r[2], par);
							}
						}
					}

					GetInfluenceArea(i_min, i_max, j_min, j_max, par.N1_v - 1, par.N2_v - 1, xs, 4, par);
					solid->Nodes[k].uf[2] = 0.0;
					for (int i = i_min; i <= i_max; ++i) {
						int i_real = i_real_v(i, par);
						for (int j = j_min; j <= j_max; ++j) {
							GeomVec r = x_v(i_real, j, par) - xs;
							solid->Nodes[k].uf[2] += v[i_real][j] * DeltaFunction(r[1], r[2], par);
							if (par.BC == periodical) {
								double x_minus = r[1] - par.L;
								double x_plus  = r[1] + par.L;
								if (x_minus < dx3) solid->Nodes[k].uf[2] += v[i_real][j] * DeltaFunction(x_minus, r[2], par);
								if (x_plus  < dx3) solid->Nodes[k].uf[2] += v[i_real][j] * DeltaFunction(x_plus , r[2], par);
							}
						}
					}

					// mass force f in Lagrange nodes
					solid->Nodes[k].f_tmp = -(solid->Nodes[k].uf - solid->Nodes[k].us) / par.d_t;
					solid->Nodes[k].f += solid->Nodes[k].f_tmp;

					double ds = solid->ds(k);
					GeomVec r = solid->Nodes[k].x / length(solid->Nodes[k].x);
					solid->Fr += dot_product(r, solid->Nodes[k].f_tmp) * ds;   // compression force applied to the solid
					solid->S += ds;
					solid->f_L += solid->Nodes[k].f_tmp * ds * par.d_x;
					r = solid->Nodes[k].x;
					r[2] += 0.;
					solid->tau_L += x_product(r, solid->Nodes[k].f_tmp) * ds * par.d_x;
				}

				std::clock_t end = std::clock();
				//std::cout << end - begin << std::endl;

				solid->Fr /= solid->S;
				for (size_t k = 0; k < solid->Nn; ++k) {
					GeomVec r = (solid->Nodes[k].x) / length(solid->Nodes[k].x);
					solid->Nodes[k].f_tmp -= solid->Fr * r;
				}

				solid->Fr_all += solid->Fr;

				CreateMatrix(Fx_temp, par.N1_u, par.N2_u);
				CreateMatrix(Fy_temp, par.N1_v, par.N2_v);

				begin = std::clock();

				for (size_t k = 0; k < solid->Nn; ++k) {

					int ix_max, ix_min;
					int jx_max, jx_min;

					int iy_max, iy_min;
					int jy_max, jy_min;

					GeomVec xs = solid->x + solid->Nodes[k].x;
					GetInfluenceArea(ix_min, ix_max, jx_min, jx_max, par.N1_u - 1, par.N2_u - 1, xs, 4, par);
					GetInfluenceArea(iy_min, iy_max, jy_min, jy_max, par.N1_v - 1, par.N2_v - 1, xs, 4, par);

					//calculating velocities us of the solid boundary
					double ds = solid->ds(k);
					double dn = sqrt(par.d_x*par.d_x * solid->Nodes[k].n[1] * solid->Nodes[k].n[1]
					               + par.d_y*par.d_y * solid->Nodes[k].n[2] * solid->Nodes[k].n[2]) / par.d_x / par.d_y;
					// calculating force force_temp for Euler nodes caused by k-th solid
					for (int i = ix_min; i <= ix_max; ++i) {
						for (int j = jx_min; j <= jx_max; ++j) {
							int i_real = i_real_u(i, par);
							GeomVec r = x_u(i_real, j, par) - xs;
							Fx_temp[i_real][j] += solid->Nodes[k].f_tmp[1] * DeltaFunction(r[1], r[2], par) * dn * ds;

							if (par.BC == periodical) {
								double x_minus = r[1] - par.L;
								double x_plus  = r[1] + par.L;
								if (x_plus  < dx3) Fx_temp[i_real][j] += solid->Nodes[k].f_tmp[1] * DeltaFunction(x_plus , r[2], par) * dn * ds;
								if (x_minus < dx3) Fx_temp[i_real][j] += solid->Nodes[k].f_tmp[1] * DeltaFunction(x_minus, r[2], par) * dn * ds;
							}
						}
					}

					for (int i = iy_min; i <= iy_max; ++i) {
						for (int j = jy_min; j <= jy_max; ++j) {
							int i_real = i_real_v(i, par);
							GeomVec r = x_v(i_real, j, par) - xs;
							Fy_temp[i_real][j] += solid->Nodes[k].f_tmp[2] * DeltaFunction(r[1], r[2], par) * dn * ds;

							if (par.BC == periodical) {
								double x_minus = r[1] - par.L;
								double x_plus  = r[1] + par.L;
								if (x_plus  < dx3) Fy_temp[i_real][j] += solid->Nodes[k].f_tmp[2] * DeltaFunction(x_plus , r[2], par) * dn * ds;
								if (x_minus < dx3) Fy_temp[i_real][j] += solid->Nodes[k].f_tmp[2] * DeltaFunction(x_minus, r[2], par) * dn * ds;
							}
						}
					}
				}

				end = std::clock();
				//std::cout << end - begin << std::endl;

				//std::cin.get();

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
						solid->f[1] += Fx_temp[i_real][j] * par.d_x * par.d_y;
						GeomVec r = x_u(i, j, par) - solid->x;
						GeomVec f;
						f[1] = Fx_temp[i_real][j];
						f[2] = 0.0;
						f[3] = 0.0;
						solid->tau += x_product(r, f) *  par.d_x * par.d_y;
					}
				}
				for (int i = iy_min; i <= iy_max; ++i) {
					for (int j = jy_min; j <= jy_max; ++j) {
						int i_real = i_real_v(i, par);
						Fy[i_real][j] += Fy_temp[i_real][j];
						solid->f[2] += Fy_temp[i_real][j] * par.d_x * par.d_y;
						GeomVec r = x_v(i, j, par) - solid->x;
						GeomVec f;
						f[1] = 0.0;
						f[2] = Fy_temp[i_real][j];
						f[3] = 0.0;
						solid->tau += x_product(r, f) *  par.d_x * par.d_y;
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
			GeomVec xs = solid.x + solid.Nodes[k].x;
			for (int i = i_min; i <= i_max; ++i) {
				for (int j = j_min; j <= j_max; ++j) {
					int i_real = i;
					if (i_real <  1      ) i_real += np1 - 2;
					if (i_real >  np1 - 2) i_real -= np1 - 2;
					GeomVec r = x_p(i_real, j, par) - xs;
					solid.Nodes[k].Eps(1, 1) += Exx[i_real][j] * DeltaFunction(r[1], r[2], par);
					solid.Nodes[k].Eps(2, 2) += Eyy[i_real][j] * DeltaFunction(r[1], r[2], par);
					solid.Nodes[k].p         +=   p[i_real][j] * DeltaFunction(r[1], r[2], par);
					if (par.BC == periodical) {
						double x_minus = r[1] - par.L;
						double x_plus  = r[1] + par.L;
						solid.Nodes[k].Eps(1, 1) += Exx[i_real][j] * DeltaFunction(x_minus, r[2], par);
						solid.Nodes[k].Eps(1, 1) += Exx[i_real][j] * DeltaFunction(x_plus , r[2], par);
						solid.Nodes[k].Eps(2, 2) += Eyy[i_real][j] * DeltaFunction(x_minus, r[2], par);
						solid.Nodes[k].Eps(2, 2) += Eyy[i_real][j] * DeltaFunction(x_plus , r[2], par);
						solid.Nodes[k].p         +=   (p[i_real][j] + dpdx_Poiseuille(par.H, par.Re)*par.L) * DeltaFunction(x_minus, r[2], par);
						solid.Nodes[k].p         +=   (p[i_real][j] - dpdx_Poiseuille(par.H, par.Re)*par.L) * DeltaFunction(x_plus , r[2], par);
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
					GeomVec r = x_c(i_real, j, par) - xs;
					solid.Nodes[k].Eps(1, 2) += Exy[i_real][j] * DeltaFunction(r[1], r[2], par);
					if (par.BC == periodical) {
						double x_minus = r[1] - par.L;
						double x_plus  = r[1] + par.L;
						solid.Nodes[k].Eps(1, 2) += Exy[i_real][j] * DeltaFunction(x_minus, r[2], par);
						solid.Nodes[k].Eps(1, 2) += Exy[i_real][j] * DeltaFunction(x_plus , r[2], par);
					}
				}
			}

		}
	}

}

void Solids_Force(std::list<Circle> &Solids, double Re) {

	for (auto& solid : Solids) {
		std::fill(solid.F_hd.begin(), solid.F_hd.end(), 0.0);
		std::fill(solid.tau_hd.begin(), solid.tau_hd.end(), 0.0);
		for (size_t k = 0; k < solid.Nn; ++k) {
			double ds = solid.ds(k);
			solid.Nodes[k].t =  ublas::prod(solid.Nodes[k].Eps, solid.Nodes[k].n) / Re    // * rho (rho = 1)
			                  - solid.Nodes[k].p  * solid.Nodes[k].n;
			solid.F_hd   += solid.Nodes[k].t * ds;
			solid.tau_hd += x_product(solid.Nodes[k].x, solid.Nodes[k].t) * ds;
		}
	}
}

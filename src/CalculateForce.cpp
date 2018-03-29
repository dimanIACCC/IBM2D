#include "stdafx.h"
#include "CalculateForce.h"
#include "Output.h"



void CalculateForce(Matrix& force_x, Matrix& force_y, std::list<Circle> &iList, Matrix& u, Matrix& v, Param par) {

	size_t nx1 = par.N1;
	size_t nx2 = par.N2 + 1;

	size_t ny1 = par.N1 + 1;
	size_t ny2 = par.N2;



	for (size_t i = 0; i < nx1; ++i) {
		for (size_t j = 0; j < nx2; ++j) {
			force_x[i][j] = 0.0;
		}
	}
	for (size_t i = 0; i < ny1; ++i) {
		for (size_t j = 0; j < ny2; ++j) {
			force_y[i][j] = 0.0;
		}
	}

	for (auto& solid : iList) {

		for (size_t k = 0; k < solid.Nn; ++k) {

			//calculating velocities us of the solid boundary
			solid.velocities();

			int i_max, i_min;
			int j_max, j_min;

			//calculating fluid velocity uf in Lagrange nodes by using near Euler nodes and discrete delta function
			GetInfluenceArea(i_min, i_max, j_min, j_max, nx1 - 1, nx1 - 1, solid.xc + solid.Nodes[k].x, 4, par);
			solid.Nodes[k].uf[1] = 0.0;
			for (int i = i_min; i <= i_max; ++i) {
				for (int j = j_min; j <= j_max; ++j) {
					int i_real = i;
					if (i_real <  0      ) i_real += nx1 - 1;
					if (i_real >  nx1 - 1) i_real -= nx1 - 1;
					if (i_real == nx1 - 1) i_real = 0;
					GeomVec xu = x_u(i_real, j, par);
					GeomVec xs = solid.xc + solid.Nodes[k].x;
					solid.Nodes[k].uf[1] += u[i_real][j] * DeltaFunction(xu[1] - xs[1], xu[2] - xs[2], par) * par.d_x * par.d_y;
					if (par.BC == periodical) {
						solid.Nodes[k].uf[1] += u[i_real][j] * DeltaFunction(xu[1] - xs[1] - par.L, xu[2] - xs[2], par) * par.d_x * par.d_y;
						solid.Nodes[k].uf[1] += u[i_real][j] * DeltaFunction(xu[1] - xs[1] + par.L, xu[2] - xs[2], par) * par.d_x * par.d_y;
					}
				}
			}

			GetInfluenceArea(i_min, i_max, j_min, j_max, ny1 - 1, ny2 - 1, solid.xc + solid.Nodes[k].x, 4, par);
			solid.Nodes[k].uf[2] = 0.0;
			for (int i = i_min; i <= i_max; ++i) {
				for (int j = j_min; j <= j_max; ++j) {
					int i_real = i;
					if (i_real <  0      ) i_real += ny1 - 2;
					if (i_real >  ny1 - 1) i_real -= ny1 - 2;
					if (i_real == 0      ) i_real  = ny1 - 2;
					if (i_real == ny1 - 1) i_real  = 1;
					GeomVec xv = x_v(i_real, j, par);
					GeomVec xs = solid.xc + solid.Nodes[k].x;
					solid.Nodes[k].uf[2] += v[i_real][j] * DeltaFunction(xv[1] - xs[1], xv[2] - xs[2], par) * par.d_x * par.d_y;
					if (par.BC == periodical) {
						solid.Nodes[k].uf[2] += v[i_real][j] * DeltaFunction(xv[1] - xs[1] - par.L, xv[2] - xs[2], par) * par.d_x * par.d_y;
						solid.Nodes[k].uf[2] += v[i_real][j] * DeltaFunction(xv[1] - xs[1] + par.L, xv[2] - xs[2], par) * par.d_x * par.d_y;
					}
				}
			}

			// mass force f in Lagrange nodes
			solid.Nodes[k].f_tmp = -(solid.Nodes[k].uf - solid.Nodes[k].us) / par.d_t;
			solid.Nodes[k].f += solid.Nodes[k].f_tmp;

			double ds = solid.ds(k);
			GeomVec r = solid.Nodes[k].x / length(solid.Nodes[k].x);
			solid.Fr += dot_product(r, solid.Nodes[k].f_tmp) * ds;
			solid.S += ds;
		}
	}

	for (auto& solid : iList) {
		solid.Fr /= solid.S;
		for (size_t k = 0; k < solid.Nn; ++k) {
			GeomVec r = (solid.Nodes[k].x) / length(solid.Nodes[k].x);
			solid.Nodes[k].f_tmp -= solid.Fr * r;
		}

		solid.Fr_all += solid.Fr;
	}



	for (auto& solid : iList) {
		CreateMatrix(force_x_temp, nx1, nx2);
		CreateMatrix(force_y_temp, ny1, ny2);
		CreateMatrix(S_x, nx1, nx2);
		CreateMatrix(S_y, ny1, ny2);

		for (size_t k = 0; k < solid.Nn; ++k) {

			int ix_max, ix_min;
			int jx_max, jx_min;

			int iy_max, iy_min;
			int jy_max, jy_min;

			GetInfluenceArea(ix_min, ix_max, jx_min, jx_max, nx1 - 1, nx2 - 1, solid.xc + solid.Nodes[k].x, 4, par);
			GetInfluenceArea(iy_min, iy_max, jy_min, jy_max, ny1 - 1, ny2 - 1, solid.xc + solid.Nodes[k].x, 4, par);

			//calculating velocities us of the solid boundary
			double ds = solid.ds(k);

			double min_w = 2.0;
			// calculating force force_temp for Euler nodes caused by k-th solid
			for (int i = ix_min; i <= ix_max; ++i) {
				for (int j = jx_min; j <= jx_max; ++j) {
					int i_real = i;
					if (i_real <  0      ) i_real += nx1 - 1;
					if (i_real >  nx1 - 1) i_real -= nx1 - 1;
					if (i_real == nx1 - 1) i_real = 0;
					
					GeomVec xu = x_u(i_real, j, par);
					GeomVec xs = solid.xc + solid.Nodes[k].x;
					double w = FunctionD((xu[1] - xs[1]) / par.d_x) * FunctionD((xu[2] - xs[2]) / par.d_y) * 4.;
					if (length(xu - xs) / par.d_x < min_w) {
						force_x_temp[i_real][j] += solid.Nodes[k].f_tmp[1] * w * ds;
						S_x[i_real][j] += ds;
					}

					if (par.BC == periodical) {

						GeomVec xu_plus = xu;
						xu_plus[1] += par.L;
						w = FunctionD((xu_plus[1] - xs[1]) / par.d_x) * FunctionD((xu_plus[2] - xs[2]) / par.d_y) * 4.;
						if (length(xu_plus - xs) / par.d_x < min_w) {
							force_x_temp[i_real][j] += solid.Nodes[k].f_tmp[1] * w * ds;
							S_x[i_real][j] += ds;
						}

						GeomVec xu_minus = xu;
						xu_minus[1] -= par.L;
						w = FunctionD((xu_minus[1] - xs[1]) / par.d_x) * FunctionD((xu_minus[2] - xs[2]) / par.d_y) * 4.;
						if (length(xu_minus - xs) / par.d_x < min_w) {
							force_x_temp[i_real][j] += solid.Nodes[k].f_tmp[1] * w * ds;
							S_x[i_real][j] += ds;
						}

					}
				}
			}

			for (int i = iy_min; i <= iy_max; ++i) {
				for (int j = jy_min; j <= jy_max; ++j) {
					int i_real = i;
					if (i_real <  0      ) i_real += ny1 - 2;
					if (i_real >  ny1 - 1) i_real -= ny1 - 2;
					if (i_real == 0      ) i_real  = ny1 - 2;
					if (i_real == ny1 - 1) i_real  = 1;
					GeomVec xv = x_v(i_real, j, par);
					GeomVec xs = solid.xc + solid.Nodes[k].x;
					double w = FunctionD((xv[1] - xs[1]) / par.d_x) * FunctionD((xv[2] - xs[2]) / par.d_y) * 4.;
					if (length(xv - xs) / par.d_y < min_w) {
						force_y_temp[i_real][j] += solid.Nodes[k].f_tmp[2] * w * ds;
						S_y[i_real][j] += ds;
					}

					if (par.BC == periodical) {

						GeomVec xv_plus = xv;
						xv_plus[1] += par.L;
						w = FunctionD((xv_plus[1] - xs[1]) / par.d_x) * FunctionD((xv_plus[2] - xs[2]) / par.d_y) * 4.;
						if (length(xv_plus - xs) / par.d_y < min_w) {
							force_y_temp[i_real][j] += solid.Nodes[k].f_tmp[2] * w * ds;
							S_y[i_real][j] += ds;
						}

						GeomVec xv_minus = xv;
						xv_minus[1] -= par.L;
						w = FunctionD((xv_minus[1] - xs[1]) / par.d_x) * FunctionD((xv_minus[2] - xs[2]) / par.d_y) * 4.;
						if (length(xv_minus - xs) / par.d_y < min_w) {
							force_y_temp[i_real][j] += solid.Nodes[k].f_tmp[2] * w * ds;
							S_y[i_real][j] += ds;
						}

					}
				}
			}
		}

		int ix_max, ix_min;
		int jx_max, jx_min;
		int iy_max, iy_min;
		int jy_max, jy_min;

		GetInfluenceArea(ix_min, ix_max, jx_min, jx_max, nx1 - 1, nx2 - 1, solid.xc, int(solid.r / par.d_x) + 4, par);
		GetInfluenceArea(iy_min, iy_max, jy_min, jy_max, ny1 - 1, ny2 - 1, solid.xc, int(solid.r / par.d_x) + 4, par);

		// summarizing force for Euler nodes and for solids
		for (int i = ix_min; i <= ix_max; ++i) {
			for (int j = jx_min; j <= jx_max; ++j) {
				int i_real = i;
				if (i_real <  0      ) i_real += nx1 - 1;
				if (i_real >  nx1 - 1) i_real -= nx1 - 1;
				if (i_real == nx1 - 1) i_real  = 0;
				if (fabs(S_x[i_real][j]) > 1.e-12) force_x_temp[i_real][j] = force_x_temp[i_real][j] / S_x[i_real][j];
				force_x[i_real][j] += force_x_temp[i_real][j];
				solid.f[1] += force_x_temp[i_real][j] * par.d_x * par.d_y;
				GeomVec r = x_u(i, j, par) - solid.xc;
				GeomVec f;
				f[1] = force_x_temp[i_real][j];
				f[2] = 0.0;
				f[3] = 0.0;
				solid.tau += x_product(r, f) *  par.d_x * par.d_y;
			}
		}
		for (int i = iy_min; i <= iy_max; ++i) {
			for (int j = jy_min; j <= jy_max; ++j) {
				int i_real = i;
				if (i_real <  0      ) i_real += ny1 - 2;
				if (i_real >  ny1 - 1) i_real -= ny1 - 2;
				if (i_real == 0      ) i_real  = ny1 - 2;
				if (i_real == ny1 - 1) i_real  = 1;
				if (fabs(S_y[i_real][j]) > 1.e-12) force_y_temp[i_real][j] = force_y_temp[i_real][j] / S_y[i_real][j];
				force_y[i_real][j] += force_y_temp[i_real][j];
				solid.f[2] += force_y_temp[i_real][j] * par.d_x * par.d_y;
				GeomVec r = x_v(i, j, par) - solid.xc;
				GeomVec f;
				f[1] = 0.0;
				f[2] = force_y_temp[i_real][j];
				f[3] = 0.0;
				solid.tau += x_product(r, f) *  par.d_x * par.d_y;
			}
		}
	}

	// copy force to non-used boundary nodes
	if (par.BC == periodical) {
		for (size_t j = 0; j < nx2; ++j) {
			force_x[nx1 - 1][j] = force_x[0][j];
		}

		for (size_t j = 0; j < ny2; ++j) {
			force_y[0][j]       = force_y[ny1 - 2][j];
			force_y[ny1 - 1][j] = force_y[1][j];
		}
	}

}

void deformation_velocity(Matrix &u, Matrix &v, Matrix &Exx, Matrix &Eyy, Matrix &Exy, Param par) {

	size_t n1 = par.N1 + 1;
	size_t n2 = par.N2 + 1;

	for (size_t i = 1; i < n1 - 1; ++i) {
		for (size_t j = 1; j < n2 - 1; ++j) {

			Exx[i][j] = (u[i][j] - u[i - 1][j]) / par.d_x;
			Eyy[i][j] = (v[i][j] - v[i][j - 1]) / par.d_y;

		}
	}

	n1 = par.N1;
	n2 = par.N2;

	for (size_t i = 0; i < n1; ++i) {
		for (size_t j = 0; j < n2; ++j) {

			Exy[i][j] = ((v[i + 1][j] - v[i][j]) / par.d_x
			           + (u[i][j + 1] - u[i][j]) / par.d_y) * 0.5;

		}
	}

}

void Solids_deformation_velocity_pressure(std::list<Circle> &Solids, Matrix &Exx, Matrix &Eyy, Matrix &Exy, Matrix &p, Param par) {

	size_t np1 = par.N1 + 1;
	size_t np2 = par.N2 + 1;

	size_t nc1 = par.N1;
	size_t nc2 = par.N2;

	for (auto& solid : Solids) {

		for (size_t k = 0; k < solid.Nn; ++k) {

			int i_max, i_min;
			int j_max, j_min;


			GetInfluenceArea(i_min, i_max, j_min, j_max, np1 - 1, np2 - 1, solid.xc + solid.Nodes[k].x, 4, par);
			solid.Nodes[k].Eps(1, 1) = 0.0;
			solid.Nodes[k].Eps(2, 2) = 0.0;
			solid.Nodes[k].p         = 0.0;
			for (int i = i_min; i <= i_max; ++i) {
				for (int j = j_min; j <= j_max; ++j) {
					int i_real = i;
					if (i_real <  0      ) i_real += np1 - 2;
					if (i_real >  np1 - 1) i_real -= np1 - 2;
					if (i_real == 0      ) i_real  = np1 - 2;
					if (i_real == np1 - 1) i_real  = 1;
					GeomVec xp = x_p(i_real, j, par);
					GeomVec xs = solid.xc + solid.Nodes[k].x;
					solid.Nodes[k].Eps(1, 1) += Exx[i_real][j] * DeltaFunction(xp[1] - xs[1], xp[2] - xs[2], par) * par.d_x * par.d_y;
					solid.Nodes[k].Eps(2, 2) += Eyy[i_real][j] * DeltaFunction(xp[1] - xs[1], xp[2] - xs[2], par) * par.d_x * par.d_y;
					solid.Nodes[k].p         +=   p[i_real][j] * DeltaFunction(xp[1] - xs[1], xp[2] - xs[2], par) * par.d_x * par.d_y;
					if (par.BC == periodical) {
						solid.Nodes[k].Eps(1, 1) += Exx[i_real][j] * DeltaFunction(xp[1] - xs[1] - par.L, xp[2] - xs[2], par) * par.d_x * par.d_y;
						solid.Nodes[k].Eps(1, 1) += Exx[i_real][j] * DeltaFunction(xp[1] - xs[1] + par.L, xp[2] - xs[2], par) * par.d_x * par.d_y;
						solid.Nodes[k].Eps(2, 2) += Eyy[i_real][j] * DeltaFunction(xp[1] - xs[1] - par.L, xp[2] - xs[2], par) * par.d_x * par.d_y;
						solid.Nodes[k].Eps(2, 2) += Eyy[i_real][j] * DeltaFunction(xp[1] - xs[1] + par.L, xp[2] - xs[2], par) * par.d_x * par.d_y;
						solid.Nodes[k].p         +=   p[i_real][j] * DeltaFunction(xp[1] - xs[1] - par.L, xp[2] - xs[2], par) * par.d_x * par.d_y;
						solid.Nodes[k].p         +=   p[i_real][j] * DeltaFunction(xp[1] - xs[1] + par.L, xp[2] - xs[2], par) * par.d_x * par.d_y;
					}
				}
			}

			GetInfluenceArea(i_min, i_max, j_min, j_max, nc1 - 1, nc2 - 1, solid.xc + solid.Nodes[k].x, 4, par);
			solid.Nodes[k].Eps(1, 2) = 0.0;
			for (int i = i_min; i <= i_max; ++i) {
				for (int j = j_min; j <= j_max; ++j) {
					int i_real = i;
					if (i_real <  0      ) i_real += nc1 - 1;
					if (i_real >  nc1 - 1) i_real -= nc1 - 1;
					if (i_real == 0      ) i_real  = nc1 - 1;
					GeomVec xc = x_c(i_real, j, par);
					GeomVec xs = solid.xc + solid.Nodes[k].x;
					solid.Nodes[k].Eps(1, 2) += Exy[i_real][j] * DeltaFunction(xc[1] - xs[1], xc[2] - xs[2], par) * par.d_x * par.d_y;
					if (par.BC == periodical) {
						solid.Nodes[k].Eps(1, 2) += Exy[i_real][j] * DeltaFunction(xc[1] - xs[1] - par.L, xc[2] - xs[2], par) * par.d_x * par.d_y;
						solid.Nodes[k].Eps(1, 2) += Exy[i_real][j] * DeltaFunction(xc[1] - xs[1] + par.L, xc[2] - xs[2], par) * par.d_x * par.d_y;
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

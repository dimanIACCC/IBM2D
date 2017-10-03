#include "CalculateForce.h"

double DeltaFunction(double x, double y, Param par) {
	return 1.0 / (par.d_x*par.d_y) * FunctionD(x / par.d_x) * FunctionD(y / par.d_y);
}

double FunctionD(double r) {
	if ((0.0 <= fabs(r)) && (fabs(r) < 1.0)) {
		return 1.0 / 8.0*(3.0 - 2.0 * fabs(r) + sqrt(1.0 + 4.0 * fabs(r) - 4.0 * r * r));
	}
	if ((1.0 <= fabs(r)) && (fabs(r) < 2.0)) {
		return 1.0 / 8.0*(5.0 - 2.0 * fabs(r) - sqrt(-7.0 + 12.0 * fabs(r) - 4.0 * r * r));
	}
	if (2.0 <= fabs(r)) {
		return 0.0;
	}
	return 0;
}

void GetInfluenceArea(size_t &i_min, size_t &i_max, size_t &j_min, size_t &j_max, size_t &Ni, size_t &Nj, GeomVec x, size_t size, Param par){
	i_max = (size_t)(x[1] / par.d_x) + size;
	i_min = (size_t)(x[1] / par.d_x) - size;

	j_max = (size_t)(x[2] / par.d_y) + size;
	j_min = (size_t)(x[2] / par.d_y) - size;

	if (i_min < 0) {
		i_min = 0;
	}
	if (j_min < 0) {
		j_min = 0;
	}
	if (i_max >= Ni) {
		i_max = Ni - 1;
	}
	if (j_max >= Nj) {
		j_max = Nj - 1;
	}
}


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
		CreateMatrix(force_x_temp, nx1, nx2);
		CreateMatrix(force_y_temp, ny1, ny2);
		for (size_t k = 0; k < solid.Nn; ++k) {

			size_t ix_max, ix_min;
			size_t jx_max, jx_min;

			size_t iy_max, iy_min;
			size_t jy_max, jy_min;

			GetInfluenceArea(ix_min, ix_max, jx_min, jx_max, nx1, nx2, solid.Nodes[k].x, 3, par);
			GetInfluenceArea(iy_min, iy_max, jy_min, jy_max, ny1, ny2, solid.Nodes[k].x, 3, par);

			//calculating velocities us of the solid boundary
			solid.velocities();

			//calculating fluid velocity uf in Lagrange nodes by using near Euler nodes and discrete delta function
			solid.Nodes[k].uf[1] = 0.0;
			for (size_t i = ix_min; i <= ix_max; ++i) {
				for (size_t j = jx_min; j <= jx_max; ++j) {
					GeomVec xu = x_u(i, j, par);
					solid.Nodes[k].uf[1] += u[i][j] * DeltaFunction(xu[1] - solid.Nodes[k].x[1], xu[2] - solid.Nodes[k].x[2], par) * par.d_x * par.d_y;
				}
			}

			solid.Nodes[k].uf[2] = 0.0;
			for (size_t i = iy_min; i <= iy_max; ++i) {
				for (size_t j = jy_min; j <= jy_max; ++j) {
					GeomVec xv = x_v(i, j, par);
					solid.Nodes[k].uf[2] += v[i][j] * DeltaFunction(xv[1] - solid.Nodes[k].x[1], xv[2] - solid.Nodes[k].x[2], par) * par.d_x * par.d_y;
				}
			}

			//calculating Integral and force f in Lagrange nodes
			solid.Nodes[k].Integral +=  (solid.Nodes[k].uf - solid.Nodes[k].us) * par.d_t;
			solid.Nodes[k].f = par.alpha_f * solid.Nodes[k].Integral + par.beta_f  *(solid.Nodes[k].uf - solid.Nodes[k].us);


			double ds = solid.ds(k);
			// calculating force force_temp for Euler nodes caused by k-th solid
			for (size_t i = ix_min; i <= ix_max; ++i) {
				for (size_t j = jx_min; j <= jx_max; ++j) {
					GeomVec xu = x_u(i, j, par);
					force_x_temp[i][j] += solid.Nodes[k].f[1] * DeltaFunction(xu[1] - solid.Nodes[k].x[1], xu[2] - solid.Nodes[k].x[2], par) * ds * ds;
				}
			}

			for (size_t i = iy_min; i <= iy_max; ++i) {
				for (size_t j = jy_min; j <= jy_max; ++j) {
					GeomVec xv = x_v(i, j, par);
					force_y_temp[i][j] += solid.Nodes[k].f[2] * DeltaFunction(xv[1] - solid.Nodes[k].x[1], xv[2] - solid.Nodes[k].x[2], par) * ds * ds;
				}
			}
		}

		size_t ix_max = 0;
		size_t ix_min = nx1;
		size_t jx_max = 0;
		size_t jx_min = nx2;

		size_t iy_max = 0;
		size_t iy_min = ny1;
		size_t jy_max = 0;
		size_t jy_min = ny2;

		for (size_t k = 0; k < solid.Nn; ++k) {

			size_t ix_max_temp, ix_min_temp;
			size_t jx_max_temp, jx_min_temp;
			size_t iy_max_temp, iy_min_temp;
			size_t jy_max_temp, jy_min_temp;

			GetInfluenceArea(ix_min_temp, ix_max_temp, jx_min_temp, jx_max_temp, nx1, nx2, solid.Nodes[k].x, 3, par);
			GetInfluenceArea(iy_min_temp, iy_max_temp, jy_min_temp, jy_max_temp, ny1, ny2, solid.Nodes[k].x, 3, par);

			if (ix_max_temp > ix_max) {
				ix_max = ix_max_temp;
			}
			if (ix_min_temp < ix_min) {
				ix_min = ix_min_temp;
			}
			if (jx_max_temp > jx_max) {
				jx_max = jx_max_temp;
			}
			if (jx_min_temp < jx_min) {
				jx_min = jx_min_temp;
			}

			if (iy_max_temp > iy_max) {
				iy_max = iy_max_temp;
			}
			if (iy_min_temp < iy_min) {
				iy_min = iy_min_temp;
			}
			if (jy_max_temp > jy_max) {
				jy_max = jy_max_temp;
			}
			if (jy_min_temp < jy_min) {
				jy_min = jy_min_temp;
			}
		}

		// summarizing force for Euler nodes and for solids
		std::fill(solid.f.begin(),   solid.f.end()  , 0.0);
		std::fill(solid.tau.begin(), solid.tau.end(), 0.0);
		for (size_t i = ix_min; i <= ix_max; ++i) {
			for (size_t j = jx_min; j <= jx_max; ++j) {
				force_x[i][j] += force_x_temp[i][j];
				solid.f[1]    += force_x_temp[i][j] * par.d_x * par.d_y;
				GeomVec r = x_u(i, j, par) - solid.xc;
				GeomVec f;
				f[1] = force_x_temp[i][j];
				f[2] = 0.0;
				f[3] = 0.0;
				solid.tau += x_product(r, f) *  par.d_x * par.d_y;
			}
		}
		for (size_t i = iy_min; i <= iy_max; ++i) {
			for (size_t j = jy_min; j <= jy_max; ++j) {
				force_y[i][j] += force_y_temp[i][j];
				solid.f[2]    += force_y_temp[i][j] * par.d_x * par.d_y;
				GeomVec r = x_v(i, j, par) - solid.xc;
				GeomVec f;
				f[1] = 0.0;
				f[2] = force_y_temp[i][j];
				f[3] = 0.0;
				solid.tau += x_product(r, f) *  par.d_x * par.d_y;
			}
		}

		if (solid.moving) {
			solid.uc    -= solid.f   * par.d_t * 1 / (solid.rho - 1) / solid.V;  // fluid density equals 1
			solid.omega -= solid.tau * par.d_t * 1 / (solid.rho - 1) / solid.I;  // angular moment I is normalized with density
		}
	}

}


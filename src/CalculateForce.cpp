#define _USE_MATH_DEFINES
#include "GeomVec.h"
#include "CalculateForce.h"

double DeltaFunction(double x, double y, Grid grid) {
	return 1.0 / (grid.d_x*grid.d_y) * FunctionD(x / grid.d_x) * FunctionD(y / grid.d_y);
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

void GetInfluenceArea(int& i_min, int& i_max, int& j_min, int& j_max, int const& Ni, int const& Nj, GeomVec x, int size, Grid grid){
	i_max = (int)((x[1] / grid.d_x) + size);
	i_min = (int)((x[1] / grid.d_x) - size);

	j_max = (int)(x[2] / grid.d_y) + size;
	j_min = (int)(x[2] / grid.d_y) - size;

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


double CalculateForce(Matrix& force_x, Matrix& force_y, list<Circle> &iList, Matrix& u, Matrix& v, Grid grid, double alpha_f, double beta_f) {

	int const nx1 = grid.N1;
	int	const nx2 = grid.N2 + 1;

	int const ny1 = grid.N1 + 1;
	int	const ny2 = grid.N2;

	for (int i = 0; i < nx1; ++i) {
		for (int j = 0; j < nx2; ++j) {
			force_x[i][j] = 0.0;
		}
	}
	for (int i = 0; i < ny1; ++i) {
		for (int j = 0; j < ny2; ++j) {
			force_y[i][j] = 0.0;
		}
	}

	for (auto& solid : iList) {
		CreateMatrix(force_x_temp, nx1, nx2);
		CreateMatrix(force_y_temp, ny1, ny2);
		for (int k = 0; k < grid.NF; ++k) {

			int ix_max, ix_min;
			int jx_max, jx_min;

			int iy_max, iy_min;
			int jy_max, jy_min;

			GetInfluenceArea(ix_min, ix_max, jx_min, jx_max, nx1, nx2, solid.Nodes[k].x, 3, grid);
			GetInfluenceArea(iy_min, iy_max, jy_min, jy_max, ny1, ny2, solid.Nodes[k].x, 3, grid);

			//calculating velocities us of the solid boundary
			solid.velocities();

			//calculating fluid velocity uf in Lagrange nodes by using near Euler nodes and discrete delta function
			solid.Nodes[k].uf[1] = 0.0;
			for (int i = ix_min; i <= ix_max; ++i) {
				for (int j = jx_min; j <= jx_max; ++j) {
					solid.Nodes[k].uf[1] += u[i][j] * DeltaFunction(i*grid.d_x - solid.Nodes[k].x[1], (j - 0.5)*grid.d_y - solid.Nodes[k].x[2], grid) * grid.d_x * grid.d_y;
				}
			}

			solid.Nodes[k].uf[2] = 0.0;
			for (int i = iy_min; i <= iy_max; ++i) {
				for (int j = jy_min; j <= jy_max; ++j) {
					solid.Nodes[k].uf[2] += v[i][j] * DeltaFunction((i - 0.5)*grid.d_x - solid.Nodes[k].x[1], j*grid.d_y - solid.Nodes[k].x[2], grid) * grid.d_x * grid.d_y;
				}
			}

			//calculating Integral and force f in Lagrange nodes
			solid.Nodes[k].Integral +=  (solid.Nodes[k].uf - solid.Nodes[k].us) * grid.d_t;
			solid.Nodes[k].f = alpha_f * solid.Nodes[k].Integral + beta_f  *(solid.Nodes[k].uf - solid.Nodes[k].us);

			// calculating force force_temp for Euler nodes caused by k-th solid
			for (int i = ix_min; i <= ix_max; ++i) {
				for (int j = jx_min; j <= jx_max; ++j) {
					force_x_temp[i][j] += solid.Nodes[k].f[1] * DeltaFunction(i*grid.d_x - solid.Nodes[k].x[1], (j - 0.5)*grid.d_y - solid.Nodes[k].x[2], grid) * solid.d_s * solid.d_s;
				}
			}

			for (int i = iy_min; i <= iy_max; ++i) {
				for (int j = jy_min; j <= jy_max; ++j) {
					force_y_temp[i][j] += solid.Nodes[k].f[2] * DeltaFunction((i - 0.5)*grid.d_x - solid.Nodes[k].x[1], j*grid.d_y - solid.Nodes[k].x[2], grid) * solid.d_s * solid.d_s;
				}
			}
		}

		int ix_max = 0;
		int ix_min = nx1;
		int jx_max = 0;
		int jx_min = nx2;

		int iy_max = 0;
		int iy_min = ny1;
		int jy_max = 0;
		int jy_min = ny2;

		for (int k = 0; k < grid.NF; ++k) {

			int ix_max_temp, ix_min_temp;
			int jx_max_temp, jx_min_temp;
			int iy_max_temp, iy_min_temp;
			int jy_max_temp, jy_min_temp;

			GetInfluenceArea(ix_min_temp, ix_max_temp, jx_min_temp, jx_max_temp, nx1, nx2, solid.Nodes[k].x, 3, grid);
			GetInfluenceArea(iy_min_temp, iy_max_temp, jy_min_temp, jy_max_temp, ny1, ny2, solid.Nodes[k].x, 3, grid);

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
		for (int i = ix_min; i <= ix_max; ++i) {
			for (int j = jx_min; j <= jx_max; ++j) {
				force_x[i][j] += force_x_temp[i][j];
				solid.f[1]    += force_x_temp[i][j] * grid.d_x * grid.d_y;
				GeomVec r, f;
				r[1] =  i        * grid.d_x - solid.xc[1];
				r[2] = (j - 0.5) * grid.d_y - solid.xc[2];
				r[3] = 0.0;
				f[1] = force_x_temp[i][j];
				f[2] = 0.0;
				f[3] = 0.0;
				solid.tau += x_product(r, f) *  grid.d_x * grid.d_y;
			}
		}
		for (int i = iy_min; i <= iy_max; ++i) {
			for (int j = jy_min; j <= jy_max; ++j) {
				force_y[i][j] += force_y_temp[i][j];
				solid.f[2]    += force_y_temp[i][j] * grid.d_x * grid.d_y;
				GeomVec r, f;
				r[1] = (i - 0.5) * grid.d_x - solid.xc[1];
				r[2] =  j        * grid.d_y - solid.xc[2];
				r[3] = 0.0;
				f[1] = 0.0;
				f[2] = force_y_temp[i][j];
				f[3] = 0.0;
				solid.tau += x_product(r, f) *  grid.d_x * grid.d_y;
			}
		}

		if (solid.moveSolid) {
			solid.uc    -= solid.f   * grid.d_t * 1 / (solid.rho - 1) / solid.V;  // fluid density equals 1
			solid.omega -= solid.tau * grid.d_t * 1 / (solid.rho - 1) / solid.I;  // angular moment I is normalized with density
		}
	}
	return 0;
}

Matrix Calculate_F_real_x(Matrix& u_n, Matrix& v_n, Matrix& u_prev, Matrix& p, Grid g, double Re) {

	double d_xx = 1.0 / (g.d_x*g.d_x);
	double d_yy = 1.0 / (g.d_y*g.d_y);
	int n1 = g.N1;
	int n2 = g.N2 + 1;
	double advective_term_n = 0.0;
	double diffusion_term = 0.0;
	double pressure = 0.0;
	double v_help = 0.0;

	CreateMatrix(result, n1, n2);

	for (int j = 1; j < (n2 - 1); ++j) {
		for (int i = 1; i < (n1 - 1); ++i) {

			v_help = 0.25 * (v_n[i][j] + v_n[i + 1][j] + v_n[i][j - 1] + v_n[i + 1][j - 1]);
			advective_term_n = u_n[i][j] * (u_n[i + 1][j] - u_n[i - 1][j]) / (2.0*g.d_x) + v_help * (u_n[i][j + 1] - u_n[i][j - 1]) / (2.0*g.d_y);
			diffusion_term = d_xx * (u_n[i + 1][j] - 2.0 * u_n[i][j] + u_n[i - 1][j]) + d_yy * (u_n[i][j + 1] - 2.0 * u_n[i][j] + u_n[i][j - 1]);
			pressure = (p[i + 1][j] - p[i][j]) / (g.d_x);
			if (j == 1 || j == n2 - 2) {
				v_help = 0.25 * (v_n[i][j] + v_n[i + 1][j] + v_n[i][j - 1] + v_n[i + 1][j - 1]);
				advective_term_n = u_n[i][j] * (u_n[i + 1][j] - u_n[i - 1][j]) / (2.0*g.d_x) + v_help * (u_n[i][j + 1] - u_n[i][j - 1]) / (1.5*g.d_y);

				if (j == 1) {
					diffusion_term = d_xx * (u_n[i + 1][j] - 2.0 * u_n[i][j] + u_n[i - 1][j]) + d_yy * (4.0*u_n[i][j + 1] - 12.0 * u_n[i][j] + 8.0*u_n[i][j - 1]) / 3.0;
				}
				if (j == n2 - 2) {
					diffusion_term = d_xx * (u_n[i + 1][j] - 2.0 * u_n[i][j] + u_n[i - 1][j]) + d_yy * (8.0*u_n[i][j + 1] - 12.0 * u_n[i][j] + 4.0*u_n[i][j - 1]) / 3.0;
				}

				result[i][j] = advective_term_n + pressure - 1.0 / (Re) * (diffusion_term);

			}

		}

	}

	// outflow du/dx = 0
	//for (int j = 0; j < n2; ++j) {

	//	// if ( j <= M2 ){
	//	// 	continue;
	//	// }

	//	//result[n1 - 1][j] = 0.0;
	//	// result[n1-1][j] = u_n[0][j];
	//	//result[0][j] = u_n[0][j];
	//	// result[0][j] = 1.0/Re * d_yy * ( u_n[0][j + 1] - 2.0 * u_n[0][j] + u_n[0][j - 1]);
	//}

	return result;


}

Matrix Calculate_F_real_y(Matrix& u_n, Matrix& v_n, Matrix& v_prev, Matrix& p, Grid g, double Re) {

	double d_xx = 1.0 / (g.d_x*g.d_x);
	double d_yy = 1.0 / (g.d_y*g.d_y);
	int n1 = g.N1 + 1;
	int n2 = g.N2;

	double advective_term_n = 0.0;
	double diffusion_term = 0.0;
	double pressure = 0.0;
	double u_help = 0.0;

	CreateMatrix(result, n1, n2);


	for (int j = 1; j < (n2 - 1); ++j) {
		for (int i = 1; i < (n1 - 1); ++i) {

			u_help           = 0.25 * (u_n[i][j] + u_n[i - 1][j] + u_n[i][j + 1] + u_n[i - 1][j + 1]);
			advective_term_n = u_help * (v_n[i + 1][j] - v_n[i - 1][j]) / (2.0*g.d_x) + v_n[i][j] * (v_n[i][j + 1] - v_n[i][j - 1]) / (2.0*g.d_y);
			diffusion_term   = d_xx * (v_n[i + 1][j] - 2.0 * v_n[i][j] + v_n[i - 1][j]) + d_yy * (v_n[i][j + 1] - 2.0 * v_n[i][j] + v_n[i][j - 1]);
			pressure         = (p[i][j + 1] - p[i][j]) / (g.d_y);

			if (i == 1 || i == n1 - 2) {
				u_help           = 0.25 * (u_n[i][j] + u_n[i - 1][j] + u_n[i][j + 1] + u_n[i - 1][j + 1]);
				advective_term_n = u_help * (v_n[i + 1][j] - v_n[i - 1][j]) / (1.5*g.d_x) + v_n[i][j] * (v_n[i][j + 1] - v_n[i][j - 1]) / (2.0*g.d_y);

				if (i == 1) {
					diffusion_term = d_xx * (4.0*v_n[i + 1][j] - 12.0 * v_n[i][j] + 8.0*v_n[i - 1][j]) / 3.0 + d_yy * (v_n[i][j + 1] - 2.0 * v_n[i][j] + v_n[i][j - 1]);
				}
				if (i == n1 - 2) {
					diffusion_term = d_xx * (8.0*v_n[i + 1][j] - 12.0 * v_n[i][j] + 4.0*v_n[i - 1][j]) / 3.0 + d_yy * (v_n[i][j + 1] - 2.0 * v_n[i][j] + v_n[i][j - 1]);
				}

				result[i][j] = advective_term_n + pressure - 1.0 / (Re) * (diffusion_term);
			}
		}

	}

	// outflow du/dx = 0
	//for (int j = 0; j < n2; ++j) {

	//	// if ( j <= M2 ){
	//	// 	continue;
	//	// }

	//	//result[n1 - 1][j] = 0.0;
	//	// result[n1-1][j] = u_n[0][j];
	//	//result[0][j] = v_n[0][j];
	//}

	return result;


}
#include "stdafx.h"
#include "PredictVel.h"

// LHS of Navier-Stokes equation operator
// LHS = (1 / d_t - 1/Re * 0.5 \Delta ) * U_new
void Calculate_A(Template &A, Param par, double Re, Direction Dir) {

	double d_u, d_v;

	if      (Dir == Du) { d_u = par.d_x;	d_v = par.d_y;	}
	else if (Dir == Dv) { d_u = par.d_y;	d_v = par.d_x;	}

	double d_uu = 1.0 / (d_u*d_u);
	double d_vv = 1.0 / (d_v*d_v);

	A.C =  1.0 /      Re *(d_uu + d_vv) + 1.0 / par.d_t;
	A.U = -1.0 / (2.0*Re)* d_vv;
	A.R = -1.0 / (2.0*Re)* d_uu;
	A.D = -1.0 / (2.0*Re)* d_vv;
	A.L = -1.0 / (2.0*Re)* d_uu;
}

// LHS of Navier-Stokes equation
Matrix Operator_Ax(Template &A, Matrix &u, Param par, Direction Dir) {

	size_t Nx = u.size();
	size_t Ny = u[0].size();

	CreateMatrix(result, Nx, Ny);

	// For direction Du  (U, R, D, L)  are  (up, right, down, left)  neighbour elements in matrix u[i][j]
	// For direction Dv  (U, R, D, L)  are  (up, right, down, left)  neighbour elements in transpose matrix (v[i][j])^T
	// In periodical problem near the boundary the neighbour element is taken from the other side of the domain
	for (size_t j = 0; j < Ny; ++j) {
		for (size_t i = 0; i < Nx; ++i) {
			result[i][j] = A.C * u[i][j]
			             + A.U * U(u, i, j, Dir, Nx, Ny)
			             + A.R * R(u, i, j, Dir, Nx, Ny)
			             + A.D * D(u, i, j, Dir, Nx, Ny)
			             + A.L * L(u, i, j, Dir, Nx, Ny);
		}
	}

	// Periodical BC
	if (par.BC == periodical) {
		if (Dir == Du) {
			for (size_t j = 1; j < Ny; ++j) {
				result[Nx - 1][j] = result[0][j];
			}
		}
		else if (Dir == Dv) {
			for (size_t j = 0; j < Ny; ++j) {
				result[0][j] = result[Nx - 2][j];
				result[Nx - 1][j] = result[1][j];
			}
		}
	}

	// Up-Down BC
	for (size_t i = 0; i < Nx; ++i) {

		// velocity
		result[i][0]      = u[i][0];
		result[i][Ny - 1] = u[i][Ny - 1];

		// derivative of the horizontal velocity
		if ((par.BC == u_infinity) && (Dir == Du)) {
			result[i][0]      = (u[i][1]      - u[i][0]     ) / par.d_y;
			result[i][Ny - 1] = (u[i][Ny - 1] - u[i][Ny - 2]) / par.d_y;
		}

	}

	// Left-Right BC
	for (size_t j = 0; j < Ny; ++j) {
		switch (par.BC) {
			case u_infinity:	result[0][j] = u[0][j];		result[Nx - 1][j] = (3.0 * u[Nx - 1][j] - 4.0 * u[Nx - 2][j] + 1.0 * u[Nx - 3][j]) / (2.0*par.d_x);		break;
			case u_inflow  :	result[0][j] = u[0][j];		result[Nx - 1][j] = (3.0 * u[Nx - 1][j] - 4.0 * u[Nx - 2][j] + 1.0 * u[Nx - 3][j]) / (2.0*par.d_x);		break;
		}
	}

	// Taylor-Green BC
	if (par.BC == Taylor_Green) {
		for (size_t i = 0; i < Nx; ++i) {
			result[i][0]      = u[i][0];
			result[i][Ny - 1] = u[i][Ny - 1];
		}
		for (size_t j = 0; j < Ny; ++j) {
			result[0][j]      = u[0][j];
			result[Nx - 1][j] = u[Nx - 1][j];
		}
	}
	
	return result;

}

// RHS of Navier-Stokes equation
Matrix CalculateB(Matrix &u_n, Matrix &v_n, Matrix &u_s, Matrix &v_s, Matrix &p, Matrix &p_new, Param par, Direction Dir) {

	size_t Nx = u_n.size();
	size_t Ny = u_n[0].size();

	double d_u, d_v;

	if      (Dir == Du) { d_u = par.d_x;	d_v = par.d_y;}
	else if (Dir == Dv) { d_u = par.d_y;	d_v = par.d_x;}

	double d_uu = 1.0 / (d_u*d_u);
	double d_vv = 1.0 / (d_v*d_v);

	double alpha = 1.;

	CreateMatrix(result, Nx, Ny);

	for (size_t i = 1; i < (Nx - 1); ++i) {
		for (size_t j = 1; j < (Ny - 1); ++j) {
			double advective_term_n    = advective_term(u_n, v_n, i, j, d_u, d_v, Dir, Nx, Ny);
			double advective_term_s    = advective_term(u_s, v_s, i, j, d_u, d_v, Dir, Nx, Ny);
			double diffusion_term_n    = diffusion_term(u_n     , i, j, d_uu, d_vv, Dir, Nx, Ny);
			double pressure_term = 0.5 * (
			                                (R(p_new, i, j, Dir, par.N1 + 1, par.N2 + 1) - p_new[i][j]) / d_u
			                              + (R(p_new, i, j, Dir, par.N1 + 1, par.N2 + 1) - p_new[i][j]) / d_u
			                             );

			result[i][j] = -(       alpha  * advective_term_n
			               + (1.0 - alpha) * advective_term_s)
			                               - pressure_term
			                               + diffusion_term_n / (2.0*par.Re)
			                               + u_n[i][j] / par.d_t;
		}
	}

	if (par.BC == periodical) {
		if (Dir == Du) {
			for (size_t j = 1; j < Ny - 1; ++j) {
				size_t i = 0;
				double advective_term_n = advective_term(u_n, v_n, i, j, d_u, d_v, Dir, par.N1, par.N2);
				double advective_term_s = advective_term(u_s, v_s, i, j, d_u, d_v, Dir, par.N1, par.N2);
				double diffusion_term_n = diffusion_term(u_n, i, j, d_uu, d_vv, Dir, par.N1, par.N2 + 1);
				double pressure_term = (0.5 * (
				                              p_new[1][j] - p_new[par.N1 - 1][j]
				                            + p_new[1][j] - p_new[par.N1 - 1][j]
				                             )
				                             - par.L * dpdx_Poiseuille(par.H, par.Re)
										) / d_u;
				result[i][j] = -(alpha * advective_term_n
				       + (1.0 - alpha) * advective_term_s)
				                       - pressure_term
				                       + diffusion_term_n / (2.0*par.Re)
				                       + u_n[i][j] / par.d_t;
				result[Nx - 1][j] = result[0][j];
			}
		}
		else if (Dir == Dv) {
			for (size_t j = 0; j < Ny - 1; ++j) {
				result[0]     [j] = result[Nx - 2][j];
				result[Nx - 1][j] = result[1][j];
			}
		}
	}


	// Up-Down BC
	for (size_t i = 0; i < Nx; ++i) {
		if      (Dir == Du) {
			result[i][0]      = 2 * par.u_wall - u_n[i][1];
			result[i][Ny - 1] = 2 * par.u_wall - u_n[i][Ny - 2];
			
			// zero derivative of the horizontal velocity
			if (par.BC == u_infinity) {
				result[i][0] = 0;
				result[i][Ny - 1] = 0;
			}
		}
		else if (Dir == Dv) {
			result[i][0] = 0;
			result[i][Ny - 1] = 0;
		}

	}

	// Left-Right BC
	for (size_t j = 0; j < Ny; ++j) {
		switch (par.BC) {
			case u_infinity:	result[0][j] = u_n[0][j]     ;		result[Nx - 1][j] = 0;		break;
			case u_inflow  :	result[0][j] = u_n[0][j]     ;		result[Nx - 1][j] = 0;		break;
		}
	}

	// Taylor-Green BC
	if (par.BC == Taylor_Green) {
		double time_exp = exp(-4 * (M_PI*M_PI) / par.Re * par.d_t * par.N_step);
		// Up-Down BC
		for (size_t i = 0; i < Nx; ++i) {
			GeomVec x_D;
			GeomVec x_U;
			if (Dir == Du) {
				x_D = x_u(i, 0     , par);
				x_U = x_u(i, Ny - 1, par);
				double u_D = sin(M_PI * x_D[1]) * cos(M_PI * x_D[2]) * time_exp;
				double u_U = sin(M_PI * x_U[1]) * cos(M_PI * x_U[2]) * time_exp;
				
				result[i][0]      = 2 * u_D - u_n[i][1];
				result[i][Ny - 1] = 2 * u_U - u_n[i][Ny - 2];
			}
			else if (Dir == Dv) {
				x_D = x_v(i, 0     , par);
				x_U = x_v(i, Ny - 1, par);
				double v_D = -sin(M_PI * x_D[2]) * cos(M_PI * x_D[1]) * time_exp;
				double v_U = -sin(M_PI * x_U[2]) * cos(M_PI * x_U[1]) * time_exp;
			}
		}
		// Left-Right BC
		for (size_t j = 0; j < Ny; ++j) {
			GeomVec x_L;
			GeomVec x_R;
			if (Dir == Du) {
				x_L = x_u(0     , j, par);
				x_R = x_u(Nx - 1, j, par);
				double u_L = sin(M_PI * x_L[1]) * cos(M_PI * x_L[2]) * time_exp;
				double u_R = sin(M_PI * x_R[1]) * cos(M_PI * x_R[2]) * time_exp;
				result[0][j]      = u_L;
				result[Nx - 1][j] = u_R;
			}
			else if (Dir == Dv) {
				x_L = x_v(0     , j, par);
				x_R = x_v(Nx - 1, j, par);
				double v_L = -sin(M_PI * x_L[2]) * cos(M_PI * x_L[1]) * time_exp;
				double v_R = -sin(M_PI * x_R[2]) * cos(M_PI * x_R[1]) * time_exp;
				result[0][j]      = 2 * v_L - u_n[1][j];
				result[Nx - 1][j] = 2 * v_R - u_n[Nx - 2][j];
			}	
		}

	}


	return result;
}

double advective_term(Matrix &u, Matrix &v, size_t i, size_t j, double d_x, double d_y, Direction Dir, size_t N1, size_t N2) {
	double v_help = 0.25 * (v[i][j] + R(v, i, j, Dir, N1, N2) + D(v, i, j, Dir, N1, N2) + RD(v, i, j, Dir, N1, N2));
	double result = u[i][j] * (R(u, i, j, Dir, N1, N2) - L(u, i, j, Dir, N1, N2)) / (2.0*d_x)
	               + v_help * (U(u, i, j, Dir, N1, N2) - D(u, i, j, Dir, N1, N2)) / (2.0*d_y);
	return result;
}

double diffusion_term(Matrix &u, size_t i, size_t j, double d_xx, double d_yy, Direction Dir, size_t N1, size_t N2) {
	double result = d_xx * (R(u, i, j, Dir, N1, N2) - 2.0 * u[i][j] + L(u, i, j, Dir, N1, N2))
	              + d_yy * (U(u, i, j, Dir, N1, N2) - 2.0 * u[i][j] + D(u, i, j, Dir, N1, N2));
	return result;
}

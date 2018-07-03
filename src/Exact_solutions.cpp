#include "Exact_solutions.h"

double ux_Poiseuille(double y, double H) {
	return (pow(H / 2.0, 2) - pow(y - H / 2.0, 2)) / pow(H / 2.0, 2);
}

double dpdx_Poiseuille(double H, double Re) {
	return 8.0 / H / H / Re;
}

double dux_dy_Poiseuille(double y, double H) {           //    dux/dy for Poiseuille flow, divided by 2, that equals angular velocity omega for non-disturbed flow
	return -(y - H / 2.0) / pow(H / 2.0, 2);
}



double exact_u(GeomVec x, Param par, double time) {
	if      (par.BC == Taylor_Green) return Taylor_Green_u(x, par.k , par.Re, time);
	else if (par.BC == Lamb_Oseen  ) return Lamb_Oseen_u  (x, par.x0, par.Re, time);
	else if (par.BC == Line_Vortex ) return Line_Vortex_u (x, par.x0, par.Re, time);
	else if (par.BC == u_infinity  ) return 1.0;
	else if (par.BC == u_inflow    ) return ux_Poiseuille(x[2], par.H);
	else if (par.BC == periodical  ) return ux_Poiseuille(x[2], par.H);
	else { std::cout << "exact_u: unknown BC" << std::endl; return 0.; }
}

double exact_v(GeomVec x, Param par, double time) {
	if      (par.BC == Taylor_Green) return Taylor_Green_v(x, par.k , par.Re, time);
	else if (par.BC == Lamb_Oseen  ) return Lamb_Oseen_v  (x, par.x0, par.Re, time);
	else if (par.BC == Line_Vortex ) return Line_Vortex_v (x, par.x0, par.Re, time);
	else if (par.BC == u_infinity  ) return 0.0;
	else if (par.BC == u_inflow    ) return 0.0;
	else if (par.BC == periodical  ) return 0.0;
	else { std::cout << "exact_v: unknown BC" << std::endl; return 0.; }
}

double exact_p(GeomVec x, Param par, double time) {
	if      (par.BC == Taylor_Green) return Taylor_Green_p(x, par.k , par.Re, time);
	else if (par.BC == Lamb_Oseen  ) return Lamb_Oseen_p  (x, par.x0, par.Re, time);
	else if (par.BC == Line_Vortex ) return Line_Vortex_p (x, par.x0, par.Re, time);
	else if (par.BC == u_infinity  ) return (par.L - x[1]) * dpdx_Poiseuille(par.H, par.Re);
	else if (par.BC == u_inflow    ) return (par.L - x[1]) * dpdx_Poiseuille(par.H, par.Re);
	else if (par.BC == periodical  ) return (par.L - x[1]) * dpdx_Poiseuille(par.H, par.Re);
	else { std::cout << "exact_p: unknown BC" << std::endl; return 0.; }
}



double Taylor_Green_u(GeomVec x, GeomVec k, double Re, double time) {
	return  sin(k[1] * x[1]) * cos(k[2] * x[2]) * Taylor_Green_time_exp(k, Re, time);
}

double Taylor_Green_v(GeomVec x, GeomVec k, double Re, double time) {
	return -sin(k[2] * x[2]) * cos(k[1] * x[1]) * Taylor_Green_time_exp(k, Re, time) * k[1] / k[2];
}

double Taylor_Green_p(GeomVec x, GeomVec k, double Re, double time) {
	double time_exp = Taylor_Green_time_exp(k, Re, time);
	return 0.5 * (pow(cos(k[2] * x[2]) * k[1] / k[2], 2) - pow(sin(k[1] * x[1]), 2)) * time_exp * time_exp;
}

double Taylor_Green_time_exp(GeomVec k, double Re, double time) {
	return exp(-(k[1] * k[1] + k[2] * k[2]) / Re * time);
}



double Lamb_Oseen_omega(double r, double Re, double time) {
	double r0 = 0.1;
	return (1 - exp(-r*r / (4 / Re * time + r0*r0))) / (2 * M_PI * r*r);
}

double Lamb_Oseen_u(GeomVec x, GeomVec x0, double Re, double time) {
	GeomVec r = x - x0;
	GeomVec omega;
	omega[3] = Lamb_Oseen_omega(length(r), Re, time);
	GeomVec uv = x_product(omega, r);
	return uv[1];
}

double Lamb_Oseen_v(GeomVec x, GeomVec x0, double Re, double time) {
	GeomVec r = x - x0;
	GeomVec omega;
	omega[3] = Lamb_Oseen_omega(length(r), Re, time);
	GeomVec uv = x_product(omega, r);
	return uv[2];
}

double Lamb_Oseen_p(GeomVec x, GeomVec x0, double Re, double time) {

	GeomVec dx = x - x0;
	double r = length(dx);
	int N = 1000;
	double dr = r / N;

	if (r < 1e-8) return 0;

	double ri = dr * 0.5;
	double p = 0;
	for (int i = 0; i <= N; i++) {
		double u = Lamb_Oseen_omega(ri, Re, time) * ri;
		p += u * u / ri * dr;
		ri += dr;
	}

	return p;
}

double Line_Vortex_u(GeomVec x, GeomVec x0, double Re, double time) {
	GeomVec dx = x - x0;
	double r = length(dx);
	double r2 = r*r;
	return   dx[2] / (2 * M_PI * r2) * exp(-0.25*Re*r2 / (time + 1));
}

double Line_Vortex_v(GeomVec x, GeomVec x0, double Re, double time) {
	GeomVec dx = x - x0;
	double r = length(dx);
	double r2 = r*r;
	return  -dx[1] / (2 * M_PI * r2) * exp(-0.25*Re*r2 / (time + 1));
}

double Line_Vortex_p(GeomVec x, GeomVec x0, double Re, double time) {
	GeomVec dx = x - x0;
	double r = length(dx);
	double r2 = r*r;
	return -Re * ei(-0.5*Re*r2 / (time + 1)) / (16 * M_PI*M_PI * (time + 1))
	           -exp(-0.5*Re*r2 / (time + 1)) / ( 8 * M_PI*M_PI * r2);
}

//function to compute p.v.int_{ -infty }^{x} exp(t)*dt / t by summation of convergent series
double ei(double x) {
	const double eu = 0.57721566490153286060651209008240243;            // Euler constant
	double term;                                                        // Term of series
	double result;

	term = x;
	result = term;
	for (int i = 2; ; i++) {
		term*= x*(i - 1) / (i*i);
		result += term;
		if (abs(term) < 1e-14) {
			break;
		}
	}
	result += eu + log(abs(x));
	return result;
}

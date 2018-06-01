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

double Taylor_Green_u(GeomVec x, double k1, double k2, double time_exp) {
	return sin(k1 * x[1]) * cos(k2 * x[2]) * time_exp;
}

double Taylor_Green_v(GeomVec x, double k1, double k2, double time_exp) {
	return -k1 / k2 * sin(k2 * x[2]) * cos(k1 * x[1]) * time_exp;
}

double Taylor_Green_p(GeomVec x, double k1, double k2, double time_exp2) {
	return 0.5 * (pow(cos(k2 * x[2]) * k1 / k2, 2) - pow(sin(k1 * x[1]), 2)) * time_exp2;
}

double Lamb_Oseen_velocity(double r, double Re, double time) {
	double r0 = 0.05;
	return (1 - exp(-r*r / (4 / Re * time + r0*r0))) / (2 * M_PI * r);
}
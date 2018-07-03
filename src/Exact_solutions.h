#pragma once
#include "GeomVec.h"
#include "Parameters.h"

double ux_Poiseuille(double y, double H);
double dux_dy_Poiseuille(double y, double H);
double dpdx_Poiseuille(double H, double Re);

double exact_u(GeomVec x, Param par, double time);
double exact_v(GeomVec x, Param par, double time);
double exact_p(GeomVec x, Param par, double time);

double Taylor_Green_u(GeomVec x, GeomVec k, double Re, double time);
double Taylor_Green_v(GeomVec x, GeomVec k, double Re, double time);
double Taylor_Green_p(GeomVec x, GeomVec k, double Re, double time);
double Taylor_Green_time_exp(GeomVec k, double Re, double time);

double Lamb_Oseen_omega(double r, double Re, double time);
double Lamb_Oseen_u(GeomVec x, GeomVec x0, double Re, double time);
double Lamb_Oseen_v(GeomVec x, GeomVec x0, double Re, double time);
double Lamb_Oseen_p(GeomVec x, GeomVec x0, double Re, double time);

double Line_Vortex_u(GeomVec x, GeomVec x0, double Re, double time);
double Line_Vortex_v(GeomVec x, GeomVec x0, double Re, double time);
double Line_Vortex_p(GeomVec x, GeomVec x0, double Re, double time);
double ei(double x);

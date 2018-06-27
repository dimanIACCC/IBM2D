#pragma once
#include "GeomVec.h"

double ux_Poiseuille(double y, double H);
double dux_dy_Poiseuille(double y, double H);
double dpdx_Poiseuille(double H, double Re);

double Taylor_Green_u(GeomVec x, double k1, double k2, double time_exp);
double Taylor_Green_v(GeomVec x, double k1, double k2, double time_exp);
double Taylor_Green_p(GeomVec x, double k1, double k2, double time_exp2);

double Lamb_Oseen_velocity(double r, double Re, double time);
double Lamb_Oseen_pressure(double r, double Re, double time);

#include "stdafx.h"
#include "Testing.h"

using namespace std;
void DoSomeTest() {
	DoTestForce();
	//
	//TO DO: realize some tests here!
}
void DoTestForce() {
	//
	//Realized test for overrunning flow on the fixed object with default input parametres
	//

	Param par;
	par.N1 *= 2.5;
	par.N2 *= 2.5;
	const double epsilon = 1e-3;

	CreateMatrix(U_n, par.N1, par.N2 + 1);
	CreateMatrix(U_new, par.N1, par.N2 + 1);
	CreateMatrix(U_prev, par.N1, par.N2 + 1);
	CreateMatrix(B_u, par.N1, par.N2 + 1);
	CreateMatrix(Force_x, par.N1, par.N2 + 1);
	CreateMatrix(V_n, par.N1 + 1, par.N2);
	CreateMatrix(V_new, par.N1 + 1, par.N2);
	CreateMatrix(V_prev, par.N1 + 1, par.N2);
	CreateMatrix(B_v, par.N1 + 1, par.N2);
	CreateMatrix(Force_y, par.N1 + 1, par.N2);
	CreateMatrix(P, par.N1 + 1, par.N2 + 1);
	CreateMatrix(Delta_P, par.N1 + 1, par.N2 + 1);
	CreateMatrix(P_Right, par.N1 + 1, par.N2 + 1);

	Matrix OperatorA_u[5];
	Matrix OperatorA_v[5];
	Calculate_A_u(OperatorA_u, par, par.Re);
	Calculate_A_v(OperatorA_v, par, par.Re);


	ApplyInitialVelocity(U_new, par); // Applying initial data to velocity 
	U_n = U_new;
	U_prev = U_new;

	std::list<Circle> solidList; // list of immersed solids
	Circle c(par.L/2, par.H/2, 0, 0, 0, par.rho, par.Nn, false, par.r);
	solidList.push_back(c);

	int n = 0; // iteration counter
	while (n <= par.N_max) {

		CalculateForce(Force_x, Force_y, solidList, U_new, V_new, par);

		//<---------- prediction of velocity --------------------------
		B_u = CalculateB_u(U_n, V_n, U_prev, V_prev, P, Force_x, par);
		B_v = CalculateB_v(U_n, V_n, U_prev, V_prev, P, Force_y, par);
#pragma omp parallel sections num_threads(2)
		{

#pragma omp section
			{
				BiCGStab(U_new, par.N1, par.N2 + 1, OperatorA_u, B_u, par, false);
			}
#pragma omp section
			{
				BiCGStab(V_new, par.N1 + 1, par.N2, OperatorA_v, B_v, par, false);
			}

		}
		//<----------end of prediction of velocity --------------------


		P_Right = Calculate_Press_Right(U_n, V_n, par);

		for (int i = 0; i < (int)Delta_P.size(); ++i) {
			for (int j = 0; j < (int)Delta_P[i].size(); ++j) {
				Delta_P[i][j] = 0.0;
			}
		}

		double eps_p = Calculate_Press_correction(Delta_P, P_Right, par, false);

		for (int i = 0; i < par.N1 + 1; ++i) {
			for (int j = 0; j < par.N2 + 1; ++j) {
				P[i][j] = P[i][j] + 0.8 * Delta_P[i][j];
			}
		}

		for (int i = 1; i < par.N1 - 1; ++i) {
			for (int j = 1; j < par.N2; ++j) {
				U_new[i][j] = U_new[i][j] - par.d_t * (Delta_P[i + 1][j] - Delta_P[i][j]) / par.d_x;
			}
		}

		for (int j = 1; j < par.N2; ++j) {
			int i = par.N1 - 1;
			U_new[i][j] = U_new[i - 1][j];
		}

		for (int i = 1; i < par.N1 + 1; ++i) {
			for (int j = 1; j < par.N2 - 1; ++j) {
				V_new[i][j] = V_new[i][j] - par.d_t * (Delta_P[i][j + 1] - Delta_P[i][j]) / par.d_y;
			}
		}

		//------------calculating eps_u--------------------------
		double eps_u = 0.0;
		for (int i = 0; i < par.N1; ++i) {
			for (int j = 0; j < par.N2 + 1; ++j) {
				if (fabs(U_n[i][j] - U_new[i][j]) > eps_u) {
					eps_u = fabs(U_n[i][j] - U_new[i][j]);
				}

				U_prev[i][j] = U_n[i][j];
				U_n[i][j] = U_new[i][j];
			}
		}
		//------------calculating eps_v--------------------------
		double eps_v = 0.0;
		for (int i = 0; i < par.N1 + 1; ++i) {
			for (int j = 0; j < par.N2; ++j) {
				if (fabs(V_n[i][j] - V_new[i][j]) > eps_v) {
					eps_v = fabs(V_n[i][j] - V_new[i][j]);
				}
				V_prev[i][j] = V_n[i][j];
				V_n[i][j] = V_new[i][j];
			}
		}
		//--------------------------------------------------------


		if (eps_u < epsilon && eps_v < epsilon) {
			double x1 = -3.15387;
			double x2 = -2.52765;
			//double x1 = -3.15;
			//double x2 = -2.53;
			int y1 = 101;
			int y2 = 202;


			if (abs((solidList.front().f[1]-x1)/(x2-x1)-(par.N1-y1)/(y2-y1))<0.5)cout << "Nice!" << abs((solidList.front().f[1] - x1) / (x2 - x1) - (par.N1 - y1) / (y2 - y1)) << endl;
			else cout << "Not good " << abs((solidList.front().f[1] - x1) / (x2 - x1) - (par.N1 - y1) / (y2 - y1)) << endl;
			cout << std::endl << "Cx = " << solidList.front().f[1] << " Cy = " << solidList.front().f[2] << std::endl;
			break;
		}


		++n;
	}

	
}
double diff(Matrix A, Matrix B) {
	double dif = 0;
	for (int i = 0; i < (int)A.size(); i++) {
		for (int j = 0; j < (int)A[0].size(); j++)
		{
			if(abs(A[i][j] - B[i][j]) > dif) dif = abs(A[i][j] - B[i][j]);
		}
	}
	return dif;
}

void ApplyInitialVelocity(Matrix &u, Param par) {

	// Poiseuille flow 
	for (int i = 0; i < par.N1; ++i) {
		for (int j = 0; j < par.N2+1; ++j) {
			u[i][j] = 1;
		}
	}
}
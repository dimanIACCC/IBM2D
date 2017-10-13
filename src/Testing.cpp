#include "stdafx.h"
#include "Testing.h"

using namespace std;
void DoTesting() {
	std::cout << "Start testing" << std::endl;
	//MakeResultDir(L"\TestsResult");
	int Re = 20;
	//DoTestForce(Re);
	Re = 100;
	DoTestForce(Re);


	std::cout << "End of testing block" << std::endl;
	
}

void DoTestForce(int Re) {
	//
	//Realized test for overrunning flow on the fixed object with default(?) input parametres and Re=20 or Re=100
	//
	bool Overflow = true;
	fs::path dir = L"TestsResult\\Overflow(Re=" + to_wstring(Re) + (wchar_t)')';
	MakeResultDir(dir);

	Param par;
	par.Re = Re;
	par.alpha_f = -8e-5;
	par.beta_f = -2e3;

	const double epsilon = 1e-7;

#pragma region BodyOfTest
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
	
	std::ofstream log;

	
	log.open(dir.append(L"\log.txt").c_str(), std::ios::out);
	SetLog(log, par);
	dir.remove_filename();
	std::ofstream forceDrug;
	std::ofstream forceLift;
	forceDrug.open(dir.append(L"forceDrug.plt").c_str(), std::ios::out);
	dir.remove_filename();
	forceLift.open(dir.append(L"forceLift.plt").c_str(), std::ios::out);
	dir.remove_filename();

	ApplyInitialVelocity(U_new, par); // Applying initial data to velocity 
	U_n = U_new;
	U_prev = U_new;

	std::list<Circle> solidList; // list of immersed solids
	//Circle c(par.L / 2, par.H / 2, 0, 0, 0, par.rho, par.Nn, false, par.r);
	Circle c(10, 3.5, 0, 0, 0, par.rho, par.Nn, false, 1);
	solidList.push_back(c);
	double minF=0, maxF=0,prevValue, curValue;
	bool increase;
	int n = 0,n0; // iteration counter
	while (n <= par.N_max) {

		CalculateForce(Force_x, Force_y, solidList, U_new, V_new, par);
		if (Re > 43) {
			CalcForceDrugLift(Force_x, n - 1, forceDrug);
			CalcForceDrugLift(Force_y, n - 1, forceLift);
		}
		n0 = 200000; //that condition depends on initial parameters and was calculated empericaly;
		if ((n > n0)&&( Re = 100)) {
			curValue = Sum(Force_x);
			if (n == n0+2) {
				(prevValue > Sum(Force_x)) ? increase = false : increase = true;
			}
			else if(n!=n0+1) {
				if ((increase == true) && (prevValue > curValue)) {
					maxF = prevValue; 
					increase = false;
				}
				if ((increase == false) && (prevValue < curValue)) {
					minF = prevValue;
					increase = true;
				}
			}
			prevValue = curValue;
			if ((maxF != 0) && (minF != 0)) {
				cout << "Amplitude is " << abs(maxF - minF) << endl;
				cout << "Average is " << abs(maxF - minF)/2 << endl;
				log << "Amplitude is " << abs(maxF - minF) << endl;
				log << "Average is " << abs(maxF - minF) / 2 << endl;
				return (void)0;
			}
		}


		//<---------- prediction of velocity --------------------------

		B_u = CalculateB_u(U_n, V_n, U_prev, V_prev, P, Force_x, par);
		B_v = CalculateB_v(U_n, V_n, U_prev, V_prev, P, Force_y, par);
#pragma omp parallel sections num_threads(2)
		{

#pragma omp section
			{
				BiCGStab(U_new, par.N1, par.N2 + 1, OperatorA_u, B_u, par, Overflow);
			}
#pragma omp section
			{
				BiCGStab(V_new, par.N1 + 1, par.N2, OperatorA_v, B_v, par, Overflow);
			}

		}
		//<----------end of prediction of velocity --------------------


		P_Right = Calculate_Press_Right(U_n, V_n, par);

		for (int i = 0; i < (int)Delta_P.size(); ++i) {
			for (int j = 0; j < (int)Delta_P[i].size(); ++j) {
				Delta_P[i][j] = 0.0;
			}
		}

		double eps_p = Calculate_Press_correction(Delta_P, P_Right, par,log, Overflow);

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

#pragma endregion

		if (eps_u < epsilon && eps_v < epsilon) {
			if (Re < 43) {
				//the line is drawn with two points (Cd1,Nx1) & (Cd2,Nx2)
				double Cd1 = -3.15387;
				double Cd2 = -2.52765;

				int Nx1 = 101;
				int Nx2 = 202;
				if (par.N1 / (double)Nx1 > 1.5) {
					if (abs((solidList.front().f[1] - Cd1) / (Cd2 - Cd1) - (par.N1 - Nx1) / (Nx2 - Nx1)) < 0.5)cout << "OK!" << endl;
					else cout << "Not OK" << endl;
					//cout << std::endl << "Cx = " << solidList.front().f[1] << " Cy = " << solidList.front().f[2] << std::endl;
				}
				else {
					double Cd_expected = Cd1 * Nx1 / (double)par.N1;
					if (abs(Cd_expected - solidList.front().f[1]) < 0.8)cout << "OK!" << endl;
					else cout << "Not OK" << endl;
				}
			}
			else cout << "That`s strange =(" << "Re =" << par.Re << endl;

			break;
		}

		if (n % par.output_step == 0) {
			std::cout << "n = " << std::setw(6) << n << "\t eps_u = " << std::fixed << eps_u << "\t eps_v = " << std::fixed << eps_v << endl;
			PushLog(log, n, eps_u, eps_v);
			log.flush();
		}

		if (n % par.output_step == 0) {
			Output(P, U_new, V_new, n, solidList, par,dir);
		}
		++n;
	}

	
}


void ApplyInitialVelocity(Matrix &u, Param par) {

	// horizantal flow 
	for (int i = 0; i < par.N1; ++i) {
		for (int j = 0; j < par.N2+1; ++j) {
			u[i][j] = 1;
		}
	}
}

double Sum(Matrix& f) {
	double sum = 0;
	for (int i = 0; i < (int)f.size(); i++) {
		for (int j = 0; j < (int)f[0].size(); j++)
		{
			sum += f[i][j];
		}
	}
	return sum;
}
void CalcForceDrugLift(Matrix& f, int n, std::ostream &stream) {
	double sum = 0;
	for (int i = 0; i < (int)f.size(); i++) {
		for (int j = 0; j < (int)f[0].size(); j++)
		{
			sum += f[i][j];
		}
	}
	stream << n << ' ' << sum << std::endl;
}

#include "stdafx.h"
#include "Testing.h"

using namespace std;
void DoTesting() {
	std::cout << "Start testing" << std::endl;
	std::cout << "=======================" << std::endl;
	Param par;
	par.BC = u_infinity;
	/*par.alpha_f = -8e-5;
	par.beta_f = -2e3;*/
	par.AddSolids_N = 0;
	fs::path dir;
	std::list<Circle> solidList; // list of immersed solids
	Circle c(10.0, 3.5, 0.0, 0.0, 0.0, par.rho, par.Nn, false, 0, 1.0);
	solidList.push_back(c);




	//1.
	//first test
	std::cout << "First test for small Re" << std::endl;
	dir = L"TestsResult\\Overflow(Re=" + to_wstring(par.Re) + (wchar_t)')' + (wchar_t)'\\';
	CreateDirectory(dir);
	par.Re = 20;
	par.WorkDir = dir.string();
	BodyOfProgram(par,solidList,true);
	std::cout << "=======================" << std::endl;

	//2.
	//second test
	std::cout << "Second test for large Re" << std::endl;
	par.Re = 100;
	dir = L"TestsResult\\Overflow(Re=" + to_wstring(par.Re) + (wchar_t)')' + (wchar_t)'\\';
	CreateDirectory(dir);
	par.WorkDir = dir.string();
	par.N_max = 200e3;	
	BodyOfProgram(par, solidList,true);


	double minF = 0, maxF = 0;
	double t1 = 0, t2 = 0;
	int n0 = 180e3;
	ifstream inForce;
	inForce.open(dir.append(L"force.plt").c_str());
	for (int i = 0; i < 3 * n0 + 5; i++) {
		string s;
		inForce >> s;
	}
	while (!inForce.eof()) {
		double number, fx,fy;
		inForce >> number >> fx >> fy;
		if (number == n0) {
			maxF = minF = fx;
			t1 = t2 = number;
		}
		if (fx > maxF) {
			maxF = fx;
			t2 = number;
		}
		if (fx < minF) {
			minF = fx;
			t1 = number;
		}
	}
	cout << "Amplitude is "; //<< abs(maxF - minF)<< endl;
	double oldAmplitude = 0.0852;
	if (abs(abs(maxF - minF) - oldAmplitude) < 0.2*oldAmplitude)	cout << "OK!" << endl;
	else cout << "Not good, difference more than 20% " << endl;
	cout << "Average is ";// << (maxF + minF) / 2 << endl;
	double oldAverage = -10.082;
	if (abs(abs(maxF + minF) / 2 - oldAverage) < 0.2*oldAverage)	cout << "OK!" << endl;
	else cout << "Not good, difference more than 20%" << endl;
	
	//не работает определение периода
	//double fs = 1 / (2 * (t2 - t1)*p.d_t); //frequency
	//cout << "Strouhal number is " << fs << endl;
	//if (abs(fs - 0.161) < 0.2*0.161)	cout << "OK!" << endl;
	//else cout << "Not good, difference more than 20%" << endl;

	std::cout << "=======================" << std::endl;


	std::cout << "End of testing block" << std::endl;
	
}


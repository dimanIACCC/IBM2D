#include "stdafx.h"
#include "Testing.h"

using namespace std;
void DoTesting() {
	std::cout << "Start testing" << std::endl;
	std::cout << "=======================" << std::endl;
	Param p;
	int Re = 20;
	fs::path dir;
	//1.
	//first test
	std::cout << "First test for small Re" << std::endl;
	dir = L"TestsResult\\Overflow(Re=" + to_wstring(Re) + (wchar_t)')' + (wchar_t)'\\';
	MakeResultDir(dir);
	BodyOfProgram(dir.string(),Re,true);
	std::cout << "=======================" << std::endl;

	//2.
	//second test
	std::cout << "Second test for large Re" << std::endl;
	Re = 100;
	dir = L"TestsResult\\Overflow(Re=" + to_wstring(Re) + (wchar_t)')' + (wchar_t)'\\';
	MakeResultDir(dir);
	BodyOfProgram(dir.string(), Re,true);

	double minF = 0, maxF = 0;
	double t1 = 0, t2 = 0;
	int n0 = 55500;//2e5;
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
	cout << "Amplitude is ";
	if (abs(abs(maxF - minF) - 0.1431) < 0.2*0.1431)	cout << "OK!" << endl;
	else cout << "Not good" << endl;
	cout << "Average is ";
	if (abs((maxF + minF) / 2 - 23.8066) < 0.2*23.8066)	cout << "OK!" << endl;
	else cout << "Not good" << endl;
	cout << "Strouhal number is ";
	double fs = 1 / (2 * (t2 - t1)*p.d_t); //frequency
	if (abs(fs - 0.161) < 0.2*0.161)	cout << "OK!" << endl;
	else cout << "Not good" << endl;

	std::cout << "=======================" << std::endl;


	std::cout << "End of testing block" << std::endl;
	
}


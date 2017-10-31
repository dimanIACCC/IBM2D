#include "stdafx.h"
#include "CalculateForce.h"
#include "Calculate_u_p.h"
#include "Output.h"


#include "Testing.h"






#pragma warning(disable : 4996)//for using <chrono>
#pragma warning(disable : 4244)//for GetInfluenceArea

// functions

<<<<<<< HEAD

void ApplyInitialData(Matrix& u, Param par);
int sgn(double x);
=======
void SetLog(std::ostream &log, Param par);
void PushLog(std::ostream &log, int n, double eps_u, double eps_v);
>>>>>>> refs/remotes/origin/Kuranakov


int main(int argc, char *argv[]) {

	std::string WorkDir = "";
	for (int i = 1; i < argc; i++) {
		std::string line, PAR, VALUE;
		line = (std::string)argv[i];
		GetParValue(line, PAR, VALUE);
		if      (PAR == "-d"){ 
		DoSomeTest();
		std::cout << "Start main program? (Y/N)" << std::endl;
		char ch;
		std::cin >> ch;
		if ((ch != 'Y') || (ch != 'y')) return 0;
		fs::path ResultFolder = L"\Result";
		MakeResultDir(ResultFolder);
		}
		else if (PAR == "-dir") if (VALUE.size() > 0) WorkDir = VALUE + '/';
	}
	
	const double epsilon = 1e-7;

	}

	Param par(WorkDir + "input.txt"); // Construct Parameters using file input.txt

	#pragma region SetMatrices
	CreateMatrix(U_n, par.N1, par.N2 + 1);
	CreateMatrix(U_new, par.N1, par.N2 + 1);
	CreateMatrix(Fx, par.N1, par.N2 + 1);
	CreateMatrix(V_n, par.N1 + 1, par.N2);
	CreateMatrix(V_new, par.N1 + 1, par.N2);
	CreateMatrix(B_v, par.N1 + 1, par.N2);
	CreateMatrix(Fy, par.N1 + 1, par.N2);
	CreateMatrix(P, par.N1 + 1, par.N2 + 1);
	ublas::matrix<Template> A_u(par.N1, par.N2 + 1);
	ublas::matrix<Template> A_v(par.N1 + 1, par.N2);
	Calculate_A(A_u, par, par.Re, Du);
	Calculate_A(A_v, par, par.Re, Dv);
	#pragma endregion SetMatrices

	std::ofstream log;
	/*std::string filelog = "Result/log.txt";
	log.open(filelog, std::ios::out);

	ApplyInitialData(U_new, par); // Applying initial data to velocity
	*/

	//Output(P, U_new, V_new, -1, solidList, par,ResultFolder);

	log.open(WorkDir + "log.txt", std::ios::out);
	SetLog(log, par);

	std::ofstream force;
	force.open(WorkDir + "force.plt", std::ios::out);
	force << "Variables = n, Fx, Fy" << std::endl;

	ApplyInitialData(U_n, P, par); // Applying initial data to velocity

	std::list<Circle> solidList; // list of immersed solids
	Read_Solids(WorkDir + "Solids.txt", solidList, par); // read Solids from file

	Output(P, U_n, V_n, -1, solidList, par, WorkDir);

	for (int n = 0; n <= par.N_max; ++n) {
		Add_Solids(solidList, n, par);

		CalculateForce(Fx, Fy, solidList, U_n, V_n, par);
		force << n << " " << Summ(Fx) << " " << Summ(Fy) << std::endl;

		Calculate_u_p(U_n, V_n, U_new, V_new, P, Fx, Fy, A_u, A_v, solidList, par, WorkDir);

		double eps_u = diff(U_n, U_new);
		double eps_v = diff(V_n, V_new);

		U_n = U_new;
		V_n = V_new;

		Solids_move(solidList, par);

		PushLog(log, n, eps_u, eps_v);
		log.flush();

		if (n % par.output_step == 0) {

			//Output(P, U_new, V_new, n, solidList, par, ResultFolder);
			Output(P, U_new, V_new, n, solidList, par, WorkDir);
		}


		const double epsilon = 1e-6;
		if (eps_u < epsilon && eps_v < epsilon) {

			//Output(P, U_new, V_new, n, solidList, par, ResultFolder);
			Output(P, U_new, V_new, n, solidList, par, WorkDir);
			break;
		}

	}
	log.close();
	std::cout << "Over" << std::endl;
	getchar();

	return 0;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



// Apply initial data for velocity
void ApplyInitialData(Matrix &u, Param par) {

	// Poiseuille flow 
	for (int i = 0; i < par.N1; ++i) {
		for (int j = 1; j < par.N2; ++j) {
			GeomVec xu = x_u(i, j, par);
			u[i][j] = ux_Poiseuille(xu[2], par.H);
		}
	}
}



int sgn(double x)
{
	(x >= 0) ? x = 1 : x = -1;
	return x;
}


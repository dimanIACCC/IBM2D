#include "CalculateForce.h"
#include "Calculate_u_p.h"
#include "Output.h"


#include "Testing.h"



#pragma warning(disable : 4996)//for using <chrono>
#pragma warning(disable : 4244)//for GetInfluenceArea

// functions

void SetLog(std::ostream &log, Param par);
void PushLog(std::ostream &log, int n, double eps_u, double eps_v);


int main(int argc, char *argv[]) {

	std::string WorkDir = "";
	for (int i = 1; i < argc; i++) {
		std::string line, PAR, VALUE;
		line = (std::string)argv[i];
		GetParValue(line, PAR, VALUE);
		if      (PAR == "-d") DoSomeTest();
		else if (PAR == "-dir") if (VALUE.size() > 0) WorkDir = VALUE + '/';
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
			Output(P, U_new, V_new, n, solidList, par, WorkDir);
		}


		const double epsilon = 1e-6;
		if (eps_u < epsilon && eps_v < epsilon) {
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
void SetLog(std::ostream& log, Param par) {

	log << "The IBM program starts.		";
	time_t t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());   // get time now
	log << ctime(&t) << std::endl;
	log << "The parameters are the following:" << std::endl;
	log << "Reynolds number               : Re  = " << par.Re << std::endl;
	log << "Channel length                : L   = " << par.L << std::endl;
	log << "Channel width                 : W   = " << par.H << std::endl;
	log << "Number of nodes on            : N1  = " << par.N1 << std::endl;
	log << "Number of nodes on            : N2  = " << par.N2 << std::endl;
	log << "Number of nodes for a particle: Nn  = " << par.Nn << std::endl;
	log << "Time step                     : tau = " << par.d_t << std::endl;
	log << "Force parameter alpha         : alpha = " << par.alpha_f << std::endl;
	log << "Force parameter beta          : beta  = " << par.beta_f << std::endl;
	log << "Tolerance for Zeidel method   : tol = " << par.Zeidel_eps << std::endl;
	log << "Step number for start of Solids adding  : AddSolids_start    = " << par.AddSolids_start    << std::endl;
	log << "Step interval for Solids adding         : AddSolids_interval = " << par.AddSolids_interval << std::endl;
	log << "Number of added Solids                  : AddSolids_N        = " << par.AddSolids_N        << std::endl;
	log << std::endl;

}

void PushLog(std::ostream& log, int n, double eps_u, double eps_v) {
	time_t t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()); // get time now
	std::string s_time = ctime(&t);
	s_time.erase(7, 1);
	s_time.erase(0, 4);
	s_time.erase(s_time.size() - 6, 5);
	log       << "n = " << std::setw(6) << n << "\t eps_u = " << std::fixed << eps_u << "\t eps_v = " << std::fixed << eps_v << "\t" << s_time;
	std::cout << "n = " << std::setw(6) << n << "\t eps_u = " << std::fixed << eps_u << "\t eps_v = " << std::fixed << eps_v << "\t" << s_time;
}

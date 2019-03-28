/*
======================================================================================================
That project is devoted to numerical simulation of disperse flows in the channels of various shapes.
======================================================================================================

#######################################################################################
# If you have some problems with launch that program, please read README.md at first. #
#######################################################################################
*/



#include "CalculateForce.h"
#include "Calculate_u_p.h"
#include "Output.h"
#include "PredictVel.h"
#include "BodyOfProgram.h"
//#include "CSR.h"


#pragma warning(disable : 4244)//for GetInfluenceArea

int main(int argc, char *argv[]) {

	//solve_pardiso();

	fs::path WorkDir = L"Result/";
	std::list<Circle> solidList; // list of immersed solids

	for (int i = 1; i < argc; i++) {
		std::string line, PAR, VALUE;
		line = (std::string)argv[i];
		GetParValue(line, PAR, VALUE);
		
		if (PAR == "-hibernation"|| PAR == "-h") {
			Matrix U_n, V_n, P;
			Param par;
			par.WorkDir = WorkDir.string();
			std::string file = par.WorkDir + VALUE;
			Awake(file, par, solidList, U_n, V_n, P);
			BodyOfProgram(par, solidList, U_n, V_n, P, false);
			return 0;
		}
		else if (PAR == "-dir") if (VALUE.size() > 0) WorkDir = VALUE + '/';

	}

	CreateDirectory(WorkDir);
	CreateDirectory(WorkDir.string() + "/Solids");


	Param par(WorkDir.string(), "input.txt");					// create the variable which contains parameters according to input data

	Read_Solids(par.WorkDir + "Solids.txt", solidList, par);	// read Solids from file and write them into list of solids

	if (par.BC == Lamb_Oseen) {
		for (auto& it : solidList) {
			it.omega_n[3] = Lamb_Oseen_omega(it.r, par.Re, 0.0    , par.Lamb_Oseen_r0);
			it.omega  [3] = Lamb_Oseen_omega(it.r, par.Re, par.d_t, par.Lamb_Oseen_r0);
		}
	}

	CreateMatrix(U_n, par.N1_u, par.N2_u);						// creation matrices for velocity
	CreateMatrix(V_n, par.N1_v, par.N2_v);						// and pressure
	CreateMatrix(P  , par.N1_p, par.N2_p);					//
	fill_exact(U_n, V_n, P, par, 0.0);                          // Initial conditions for velocity and pressure

	BodyOfProgram(par, solidList, U_n, V_n, P);					// start solver


	std::cout << "The End" << std::endl;
	getchar();

	return 0;
}


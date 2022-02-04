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
			BodyOfProgram(par, solidList, U_n, V_n, P);
			return 0;
		}
		else if (PAR == "-post") {
			Param par(VALUE, "/input.txt");					// create the variable which contains parameters according to input data
			par.WorkDir = VALUE + "/";
			std::string file = "history_post_average10_dist5";
			history_init(par.WorkDir, file, par.BC);

			int i = 1;
			while (Read_plt(par.WorkDir + "step" + std::to_string(i * 100) + ".plt", par, solidList)) {
				std::cout << "i = " << i*100 << ", start new file" << std::endl;

				double h_average;
				h_average_of_Solids_Layer(solidList, par, h_average);
				history_log(par.WorkDir, file, par.N_step*par.d_t, h_average, 0., 0.);
				solidList.clear();
				i++;
			}
			std::cin.ignore();

			return 0;
		}
		else if (PAR == "-dir") if (VALUE.size() > 0) WorkDir = VALUE + '/';

	}

	CreateDir(WorkDir);
	CreateDir(WorkDir.string() + "/Solids");

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
	CreateMatrix(P  , par.N1_p, par.N2_p);					    //
	fill_exact(U_n, V_n, P, par, 0.0, par.d_t*0.5);             // Initial conditions for velocity and pressure

	BodyOfProgram(par, solidList, U_n, V_n, P);					// start solver


	std::cout << "The End" << std::endl;
	getchar();

	return 0;
}


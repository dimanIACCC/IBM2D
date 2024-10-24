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

#pragma warning(disable : 4244)//for GetInfluenceArea

int main(int argc, char *argv[]) {

	fs::path WorkDir = L"Result";
	std::string file = "input.txt";

	Param par;
	std::vector<Solid> Solids; // list of immersed solids
	std::vector<Node> Nodes; // list of Nodes
	Matrix U_n, V_n, P_n;

	for (int i = 1; i < argc; i++) {
		std::string line, PAR, VALUE;
		line = (std::string)argv[i];
		GetParValue(line, PAR, VALUE);

		if (PAR == "-dir") {
			if (VALUE.size() > 0) WorkDir = VALUE;
		}
		else if (PAR == "-h") {
			if (VALUE.size() > 0) file = VALUE;
		}
		else if (PAR == "-post") {
			if (VALUE.size() > 0) WorkDir = VALUE;
			Post(WorkDir, par, Solids, Nodes, U_n, V_n, P_n);
			return 0;
		}
	}

	CreateDir(WorkDir);
	CreateDir(WorkDir.string() + "/" + "Solids");
	file = WorkDir.string() + "/" + file;
	par.WorkDir = WorkDir.string();

	Load_Data(file, par, Solids, Nodes, U_n, V_n, P_n);

	BodyOfProgram(par, Solids, Nodes, U_n, V_n, P_n);

	return 0;
}


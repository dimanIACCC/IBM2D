#include "stdafx.h"
#include "CalculateForce.h"
#include "Calculate_u_p.h"
#include "Output.h"
#include "PredictVel.h"
#include "Testing.h"
#include "BodyOfProgram.h"





#pragma warning(disable : 4996)//for using <chrono>
#pragma warning(disable : 4244)//for GetInfluenceArea



int main(int argc, char *argv[]) {

	fs::path WorkDir = L"\Result\\";
	for (int i = 1; i < argc; i++) {
		std::string line, PAR, VALUE;
		line = (std::string)argv[i];
		GetParValue(line, PAR, VALUE);
		if      (PAR == "-d"){ 
		 DoTesting();
		 std::cout << "Start main program? (Y/N)" << std::endl;
		 char ch = std::cin.get();
		 if ((ch != 'Y') || (ch != 'y')) return 0;

		}
		else if (PAR == "-dir") if (VALUE.size() > 0) WorkDir = VALUE + '/';
	}

	CreateDirectory(WorkDir);

	Param par(WorkDir.string(), "input.txt");
	std::list<Circle> solidList; // list of immersed solids
	Read_Solids(par.WorkDir + "Solids.txt", solidList, par); // read Solids from file

	BodyOfProgram(par, solidList);
	

	std::cout << "The End" << std::endl;
	getchar();

	return 0;
}



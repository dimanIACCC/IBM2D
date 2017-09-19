#include "stdafx.h"
#include "CalculateForce.h"
#include "Output.h"
#include "BiCGStab.h"
#include "Calculate_press.h"
#include "PredictVel.h"
//#include "Testing.h"



#pragma warning(disable : 4996)//for using <chrono>
#pragma warning(disable : 4244)//for GetInfluenceArea

// fuctions

void SetLog(std::ostream &log, Param par);
void PushLog(std::ostream &log, int n, double eps_u, double eps_v);
void ApplyInitialData(Matrix& u, Param par);

int sgn(double x);


int main(int argc, char *argv[]) {
/*
	for (int i = 1; i < argc; i++) {
		if((std::string)argv[i] == (std::string)"-d") DoSomeTest();
}*/

	const double epsilon = 1e-3;

	Param par("input.txt"); // Construct Parameters using file input.txt
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
	//-----------creating Result folder --------------
	//char current_work_dir[FILENAME_MAX];
	//_getcwd(current_work_dir, sizeof(current_work_dir));
	//strcat_s(current_work_dir, "\\Result");
	//_mkdir(current_work_dir);
	//-------------------------------------------------
	std::string filelog = "Result/log.txt";
	log.open(filelog, std::ios::out);
	SetLog(log, par);

	ApplyInitialData(U_new, par); // Applying initial data to velocity 
	U_n = U_new;
	U_prev = U_new;

	std::list<Circle> solidList; // list of immersed solids
	Read_Solids("Solids.txt", solidList, par); // read Solids from file

	Output(P, U_new, V_new, -1, solidList, par);

	int n = 0; // iteration counter
	while (n <= par.N_max) {

		Add_Solids(solidList, 20, n, 0, 200, par);

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
		//ExplicPredVel(U_new,V_new,U_n,V_n,P,Force_x,Force_y,par);

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


		//--------------COLLISION CHECK---------------------------
		///--------------collisions between particles---------------------------------
		for (auto one = solidList.begin(); one != solidList.end(); one++) {
			for (auto two = next(one); two != solidList.end(); two++) {
				if (Collide(*one, *two, par)) if (Debug) std::cout << "Collision detected" << std::endl;
			}
		}
		///--------------end of collisions between particles---------------------------------

		///-------------collision with walls--------------------------
		for (auto one = solidList.begin(); one != solidList.end(); one++) {
			double DistUpper = par.H - one->xc[2];//<----distance to upper wall
			double DistLower = one->xc[2];//<-------distance to lower wall
			if (DistUpper < par.k_dist * one->r || DistLower < par.k_dist * one->r) {
				if (Debug) std::cout << "Collision with wall detected" << std::endl;
				one->uc[2] = -one->uc[2];
			}
		}
		///-------------end of collision with walls--------------------------

		//--------------END OF COLLISION CHECK---------------------------------------------------


		for (auto it = solidList.begin(); it != solidList.end();) {

			it->move(par.d_t);

			//delete bodies which move 95% of length
			if (it->xc[1] > par.L*0.95) {
				solidList.erase(it++);
			}
			else {
				++it;
			}
		}

		PushLog(log, n, eps_u, eps_v);
		log.flush();

		if (n % par.output_step == 0) {
			Output(P, U_new, V_new, n, solidList, par);
		}



		if (eps_u < epsilon && eps_v < epsilon) {
			Output(P, U_new, V_new, n, solidList, par);
			break;
		}


		++n;
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
	log << std::endl;

}

void PushLog(std::ostream& log, int n, double eps_u, double eps_v) {
	time_t t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()); // get time now
	std::string s_time = ctime(&t);
	s_time.erase(7, 1);
	s_time.erase(0, 4);
	s_time.erase(s_time.size() - 6, 5);
	log << "n = " << std::setw(6) << n << "\t eps_u = " << std::fixed << eps_u << "\t eps_v = " << std::fixed << eps_v << "\t" << s_time;
	std::cout << "n = " << std::setw(6) << n << "\t eps_u = " << std::fixed << eps_u << "\t eps_v = " << std::fixed << eps_v << "\t" << s_time;
}

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

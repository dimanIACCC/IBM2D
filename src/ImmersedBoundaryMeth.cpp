#include "stdafx.h"
#include "SolidBody.h"
#include "Grid.h"
#include "Calculate_A.h"
#include "CalculateForce.h"
#include "Calculate_B.h"
#include "Output.h"
#include "BiCGStab.h"
#include "Calculate_press.h"


#pragma warning(disable : 4996)//for using <chrono>
#pragma warning(disable : 4244)//for GetInfluenceArea

#define CreateMatrix(name, n, m) Matrix name(n,std::vector<double>(m, 0))

bool Debug = false;
bool InelasticCollision = false; //Perfectly inelastic collision --- абсолютно неупругие столкновения

using namespace std;


// fuctions

void InputData(Grid& grid, double &M, int &Re, double &alpha_f, double &beta_f, double& Zeidel_eps, int& output_step, int& N_max, int& N_Zeidel);
void SetLog(ostream &log, Grid grid, double M, double Re, double alpha_f, double beta_f, double Zeidel_eps);
void PushLog(ostream &log, int n, double eps_u, double eps_v);
void ApplyInitialData(Matrix& u, Grid grid);
int sgn(double x);



int main() {
	int Re;
	int N_max = 0; // number of total iterations
	double alpha_f;
	double beta_f;
	int N_Zeidel; // Number of iterations in Zeidel method
	double Zeidel_eps;
	double Cd; //drag coefficent
	double Cl; // lift coefficent
	double x = 0.0;
	double y = 0.0;
	double r = 0.5;
	double m;
	const double epsilon = 1e-3;
	int output_step = 0; //frequency of output

	// declaring variables
	Grid grid;
	double eps_u = 0.0;
	double eps_v = 0.0;
	double eps_p = 0.0;



	int n = 0; // iteration counter
	InputData(grid, m, Re, alpha_f, beta_f, Zeidel_eps, output_step, N_max, N_Zeidel); // Get value of some variables
	CreateMatrix(U_n, grid.N1, grid.N2 + 1);
	CreateMatrix(U_new, grid.N1, grid.N2 + 1);
	CreateMatrix(U_prev, grid.N1, grid.N2 + 1);
	CreateMatrix(B_u, grid.N1, grid.N2 + 1);
	CreateMatrix(Force_x, grid.N1, grid.N2 + 1);
	CreateMatrix(V_n, grid.N1 + 1, grid.N2);
	CreateMatrix(V_new, grid.N1 + 1, grid.N2);
	CreateMatrix(V_prev, grid.N1 + 1, grid.N2);
	CreateMatrix(B_v, grid.N1 + 1, grid.N2);
	CreateMatrix(Force_y, grid.N1 + 1, grid.N2);
	CreateMatrix(P, grid.N1 + 1, grid.N2 + 1);
	CreateMatrix(Delta_P, grid.N1 + 1, grid.N2 + 1);
	CreateMatrix(P_Right, grid.N1 + 1, grid.N2 + 1);

	Matrix OperatorA_u[5];
	for (int i = 0; i < 5; i++) {
		OperatorA_u[i].resize(grid.N1);
		for (int j = 0; j < grid.N1; j++) {
			OperatorA_u[i][j].resize(grid.N2 + 1);
			fill(OperatorA_u[i][j].begin(), OperatorA_u[i][j].end(), 0);
		}
	}

	Matrix OperatorA_v[5];
	for (int i = 0; i < 5; i++) {
		OperatorA_v[i].resize(grid.N1 + 1);
		for (int j = 0; j < grid.N1 + 1; j++) {
			OperatorA_v[i][j].resize(grid.N2);
			fill(OperatorA_v[i][j].begin(), OperatorA_v[i][j].end(), 0);
		}
	}
	Calculate_A_u(OperatorA_u, grid, Re);
	Calculate_A_v(OperatorA_v, grid, Re);


	// list of immersed solids
	list<Circle> solidList;
	ofstream output; // for Drag and Lift coefficents
	ofstream press_output; // press
	ofstream log;
	string filename = "Result/coefficent.plt";
	//string filepress = "Result/eps_pressure.plt";
	string filelog = "Result/log.txt";
	log.open(filelog, ios::out);
	SetLog(log, grid, m, Re, alpha_f, beta_f, Zeidel_eps);
	log << endl;

	
	ApplyInitialData(U_new, grid); // Applying initial data to velocity 
	U_n = U_new;
	U_prev = U_new;


	//Firstly adding some circles
	if (!Debug) {
		Circle c1(3.5, 2.1, r, n, grid);
		Circle c2(3.5, 4.9, r, n, grid);
		Circle c3(1.5, 1.9, r, n, grid);
		Circle c4(1.5, 5.1, r, n, grid);
		solidList.push_back(c1);
		solidList.push_back(c2);
		solidList.push_back(c3);
		solidList.push_back(c4);
	}
	else {
		/*Circle c1(1, 5, r, n, grid);
		solidList.push_back(c1);
		solidList.begin()->U = 1;
		solidList.begin()->V = -1;
		Circle c2(1, 1, r, n, grid);
		solidList.push_back(c2);
		next(solidList.begin())->U = 1;
		next(solidList.begin())->V = 1;
		eps_u = 1.0;
		eps_v = 1.0;*/
	}



	CalculateForce_X(Force_x, solidList, U_new, r, Cd, grid, alpha_f, beta_f, m);
	CalculateForce_Y(Force_y, solidList, V_new, r, Cl, grid, alpha_f, beta_f, m);

	OutputVelocity_U(U_new, -1, output_step, solidList, grid);
	OutputVelocity_V(V_new, -1, output_step, solidList, grid);

	while (n <= N_max) {

		//creation new solids
		if (!Debug) {
			if (n > 0 && fmod(n*grid.d_t, 1.5) == 0.0) {
				//double tmp = fmod(1200*grid.d_t, 1.5);
				//if (n > 0){
				int chance = 80;
				int rnd;
				rnd = rand() % 100 + 1;
				if (rnd <= chance) {
					x = 1 + ((rand() % 100 + 1) / 100.0);
					y = 1 + ((rand() % 200 + 1) / 100.0);
					Circle c(x, y, r, n, grid);
					solidList.push_back(c);
				}
				rnd = rand() % 100 + 1;
				if (rnd <= chance) {
					x = 1 + ((rand() % 100 + 1) / 100.0);
					y = 4 + ((rand() % 200 + 1) / 100.0);
					Circle c(x, y, r, n, grid);
					solidList.push_back(c);
				}

				rnd = rand() % 100 + 1;
				if (rnd <= chance) {
					x = 3 + ((rand() % 100 + 1) / 100.0);
					y = 1 + ((rand() % 200 + 1) / 100.0);
					Circle c(x, y, r, n, grid);
					solidList.push_back(c);
				}

				rnd = rand() % 100 + 1;
				if (rnd <= chance) {
					x = 3 + ((rand() % 100 + 1) / 100.0);
					y = 4 + ((rand() % 200 + 1) / 100.0);
					Circle c(x, y, r, n, grid);
					solidList.push_back(c);
				}
			}
			}
			eps_u = 0.0;
			eps_v = 0.0;

			//<---------- prediction of velocity --------------------------
			B_u = CalculateB_u(U_n, V_n, U_prev, V_prev, P, Force_x, grid, Re);
			B_v = CalculateB_v(U_n, V_n, U_prev, V_prev, P, Force_y, grid, Re);

			BiCGStab(U_new, grid.N1, grid.N2 + 1, OperatorA_u, B_u, grid);
			BiCGStab(V_new, grid.N1 + 1, grid.N2, OperatorA_v, B_v, grid);
			//<----------end of prediction of velocity --------------------------


			P_Right = Calculate_Press_Right(U_n, V_n, grid);

			for (int i = 0; i < (int)Delta_P.size(); ++i)
				for (int j = 0; j < (int)Delta_P[i].size(); ++j) {
					Delta_P[i][j] = 0.0;
				}


			eps_p = Calculate_Press_correction(Delta_P, P_Right, N_Zeidel, Zeidel_eps, grid);
			press_output << n << ' ' << eps_p << endl; //<---writing in closed stream



			for (int i = 0; i < grid.N1 + 1; ++i) {
				for (int j = 0; j < grid.N2 + 1; ++j) {
					P[i][j] = P[i][j] + 1.0 * Delta_P[i][j]; // 0.8 changed to 1.0
				}
			}

			for (int i = 1; i < grid.N1 - 1; ++i) {
				for (int j = 1; j < grid.N2; ++j) {
					U_new[i][j] = U_new[i][j] - grid.d_t * (Delta_P[i + 1][j] - Delta_P[i][j]) / grid.d_x;
				}
			}

			for (int j = 1; j < grid.N2; ++j) {
				int i = grid.N1 - 1;
				U_new[i][j] = U_new[i - 1][j];
			}

			for (int i = 1; i < grid.N1 + 1; ++i) {
				for (int j = 1; j < grid.N2 - 1; ++j) {
					V_new[i][j] = V_new[i][j] - grid.d_t * (Delta_P[i][j + 1] - Delta_P[i][j]) / grid.d_y;
				}
			}

			//calculaating eps_u
			for (int i = 0; i < grid.N1; ++i) {
				for (int j = 0; j < grid.N2 + 1; ++j) {
					if (fabs(U_n[i][j] - U_new[i][j]) > eps_u) {
						eps_u = fabs(U_n[i][j] - U_new[i][j]);
					}

					U_prev[i][j] = U_n[i][j];
					U_n[i][j] = U_new[i][j];
				}
			}
			//calculaating eps_v
			for (int i = 0; i < grid.N1 + 1; ++i) {
				for (int j = 0; j < grid.N2; ++j) {
					if (fabs(V_n[i][j] - V_new[i][j]) > eps_v) {
						eps_v = fabs(V_n[i][j] - V_new[i][j]);
					}
					V_prev[i][j] = V_n[i][j];
					V_n[i][j] = V_new[i][j];
				}
			}


			CalculateForce_X(Force_x, solidList, U_new, r, Cd, grid, alpha_f, beta_f, m);
			CalculateForce_Y(Force_y, solidList, V_new, r, Cl, grid, alpha_f, beta_f, m);
		


//--------------COLLISION CHECK---------------------------
		list<Circle>::iterator first;
		list<Circle>::iterator second;
		///--------------collisions between particles---------------------------------
		for (auto one = solidList.begin(); one != solidList.end(); one++) {
			for (auto two = next(one); two != solidList.end(); two++) {
				double distance = sqrt(pow(one->x - two->x, 2) + pow(one->y - two->y, 2));//<----distance between two particles 
				if (one->x < two->x) {
					first = one;
					second = two;
				}
				else
				{
					first = two;
					second = one;
				}
				if (distance <= (2 * first->r)) {
					if (Debug) cout << "COLLISION DETECTED";
					if (InelasticCollision) {
						//Perfectly inelastic collision
						first->U = (first->U + second->U) / 2;
						first->V = (first->V + second->V) / 2;
						second->U = first->U;
						second->V = first->V;

					}
					else {
						//this is perfectly elastic impact
						double v1, v2;
						v1 = sgn(first->U)*sqrt(pow(first->V, 2) + pow(first->U, 2));
						v2 = sgn(second->U)*sqrt(pow(second->V, 2) + pow(second->U, 2));

						double tetta_1, tetta_2, phi;
						if (v1 != 0) { tetta_1 = atan(first->V / first->U); }
						else { tetta_1 = 0; } //??if (v1 == 0) { tetta_1 = M_PI_2 - tetta_2 / 2; }
						if (v2 != 0) { tetta_2 = atan(second->V / second->U); }
						else { tetta_2 = 0; } //??if (v2 == 0) { tetta_2 = M_PI_2 - tetta_1 / 2; }
						phi = atan((second->y - first->y) / (second->x - first->x));
						first->U = v2*cos(tetta_2 - phi)*cos(phi) + v1*sin(tetta_1 - phi)*cos(phi + M_PI_2);
						first->V = v2*cos(tetta_2 - phi)*sin(phi) + v1*sin(tetta_1 - phi)*sin(phi + M_PI_2);
						second->U = v1*cos(tetta_1 - phi)*cos(phi) + v2*sin(tetta_2 - phi)*cos(phi + M_PI_2);
						second->V = v1*cos(tetta_1 - phi)*sin(phi) + v2*sin(tetta_2 - phi)*sin(phi + M_PI_2);
					}
				}
			}
		}
		///--------------end of collisions between particles---------------------------------

		///-------------collision with walls--------------------------
		for (auto one = solidList.begin(); one != solidList.end(); one++) {
			double DistUpper = grid.H - one->y;//<----distance to upper wall
			double DistLower = one->y;//<-------distance to lower wall
			if (DistUpper < one->r || DistLower < one->r) {
				if (Debug) cout << "COLLISION DETECTED";
				one->V = -one->V;
			}
		}
		///-------------end of collision with walls--------------------------

//--------------END OF COLLISION CHECK---------------------------------------------------


		for (auto it = solidList.begin(); it != solidList.end();) {
			if ((it->moveSolid == false)) {// && (n - it->start_n > 10)) { 
				it->moveSolid = true;
			}

			if (it->moveSolid) {
				//update position
				for (int k = 0; k < grid.NF; ++k) {
					it->Bound[0][k] += it->U * grid.d_t;
					it->Bound[1][k] += it->V * grid.d_t;
				}
				it->x += it->U * grid.d_t;
				it->y += it->V * grid.d_t;
			}

			//delete bodies which move 95% of length
			if (it->x > grid.L*0.95) {
				solidList.erase(it++);
			}
			else {
				++it;
			}
		}



		if (n < 1000 || (n > 1000 && 0 == n % 100)) {
			std::cout << "n  = " << n << " | eps_u = " << eps_u << " | eps_v = " << eps_v << "		";
			time_t t = chrono::system_clock::to_time_t(chrono::system_clock::now());   // get time now
			string s_time = ctime(&t);
			s_time.erase(7, 1);
			s_time.erase(0, 4);
			s_time.erase(s_time.size() - 6, 5);
			std::cout << s_time;

			PushLog(log, n, eps_u, eps_v);
			log.flush();
		}

		if (0 == n % output_step) {
			OutputVelocity_U(U_new, n, output_step, solidList, grid);
			OutputVelocity_V(V_new, n, output_step, solidList, grid);
			OutputPressure(P, n, output_step, solidList, grid);
		}



		if (eps_u < epsilon && eps_v < epsilon) {
			OutputVelocity_U(U_new, n, output_step, solidList, grid);
			OutputVelocity_V(V_new, n, output_step, solidList, grid);
			OutputPressure(P, n, output_step, solidList, grid);
			break;
		}


		++n;
	}
	log.close();
	std::cout << "Over" << endl;
	getchar();

	return 0;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SetLog(ostream& log, Grid grid, double M, double Re, double alpha_f, double beta_f, double Zeidel_eps) {

	log << "The IBM program starts.		";
	time_t t = chrono::system_clock::to_time_t(chrono::system_clock::now());   // get time now
	log << ctime(&t) << endl;
	log << "The parameters are the following:" << endl;
	log << "Mass of a particle            : M   = " << M << endl;
	log << "Reynolds number               : Re  = " << Re << endl;
	log << "Channel length                : L   = " << grid.L << endl;
	log << "Channel width                 : W   = " << grid.H << endl;
	log << "Number of nodes on            : Nх  = " << grid.N1 << endl;
	log << "Number of nodes on            : Ny  = " << grid.N2 << endl;
	log << "Number of nodes for a particle: Np  = " << grid.NF << endl;
	log << "Time step                     : tau = " << grid.d_t << endl;
	log << "Force parameter alpha         : alpha = " << alpha_f << endl;
	log << "Force parameter beta          : beta  = " << beta_f << endl;
	log << "Tolerance for Zeidel method   : tol = " << Zeidel_eps << endl;

}

void PushLog(ostream& log, int n, double eps_u, double eps_v) {
	time_t t = chrono::system_clock::to_time_t(chrono::system_clock::now());
	log << "n  = " << n << " | eps_u = " << eps_u << " | eps_v = " << eps_v << '\t';
	string s_time = ctime(&t);
	s_time.erase(7, 1);
	s_time.erase(0, 4);
	s_time.erase(s_time.size() - 6, 5);
	log << s_time;
}

// Apply initial data for velocity
void ApplyInitialData(Matrix &u, Grid grid) {

	// Poiseuille flow 
	for (int i = 0; i < grid.N1; ++i) {
		for (int j = 1; j < grid.N2; ++j) {
			u[i][j] = (pow((grid.H) / 2.0, 2) - pow((j - 0.5)*grid.d_y - grid.H / 2.0, 2));
		}
	}
}

int sgn(double x)
{
	(x >= 0) ? x = 1 : x = -1;
	return x;
}

void InputData(Grid& grid, double &M, int &Re, double &alpha_f, double &beta_f, double& Zeidel_eps, int& output_step, int& N_max, int& N_Zeidel) {

	ifstream input;
	string filename = "input.txt";
	input.open(filename.c_str());

	input >> M;
	input >> Re;
	input >> grid.L;
	input >> grid.H;
	input >> grid.N1;
	input >> grid.N2;
	input >> grid.NF;
	input >> grid.d_t;
	input >> alpha_f;
	input >> beta_f;
	input >> output_step;
	input >> N_max;
	input >> N_Zeidel;
	input >> Zeidel_eps;
	input.close();

	grid.d_x = grid.L / (grid.N1 - 1);
	grid.d_y = grid.H / (grid.N2 - 1);

}

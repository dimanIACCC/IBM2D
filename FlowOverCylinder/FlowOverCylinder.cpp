// FlowOverCylinder.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include "SolidBody.h"
#include "Grid.h"
#include "CalculateForce.h"
#include "Output.h"
#include "BiCGStab.h"
#include "Calculate_press.h"
#include "PredictVel.h"




#pragma warning(disable : 4996)//for using <chrono>
#pragma warning(disable : 4244)//for GetInfluenceArea

#define CreateMatrix(name, n, m) Matrix name(n,std::vector<double>(m, 0))

bool InelasticCollision = false; //Perfectly inelastic collision --- абсолютно неупругие столкновения

using namespace std;



// fuctions

Circle InputData(Grid& grid, double &M, int &Re, double &alpha_f, double &beta_f, double& Zeidel_eps, int& output_step, int& N_max, int& N_Zeidel);
double Summ(Matrix& force);
void SetLog(ostream &log, Grid grid, double M, double Re, double alpha_f, double beta_f, double Zeidel_eps);
void PushLog(ostream &log, int n, double eps_u, double eps_v);
void ApplyInitialData(Matrix& u, Grid grid);
int sgn(double x);



int main() {
#pragma region variables
	int Re;
	int N_max = 0; // number of total iterations
	double alpha_f;
	double beta_f;
	int N_Zeidel; // Number of iterations in Zeidel method
	double Zeidel_eps;
	double Cd; //drag coefficent
	double Cl; // lift coefficent
	double r = 0.5;
	double m;
	const double epsilon = 1e-7;
	int output_step = 0; //frequency of output
	Grid grid;
	double eps_u = 0.0;
	double eps_v = 0.0;
	double eps_p = 0.0;
	list<Circle> solidList;
#pragma endregion variables
	Circle cylinder(0, 0, 0, 0, grid);
	int n = 0; // iteration counter
	cylinder = InputData(grid, m, Re, alpha_f, beta_f, Zeidel_eps, output_step, N_max, N_Zeidel); // Get value of some variables
	solidList.push_back(cylinder);

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

	ofstream output; // for Drag and Lift coefficents
	ofstream press_output; // press
	ofstream log;
	ofstream fDrug, fLift;
	//-----------creating Result folder --------------
	/*char current_work_dir[FILENAME_MAX];
	_getcwd(current_work_dir, sizeof(current_work_dir)); //not crossplatform solution
	strcat_s(current_work_dir, "\\Result");
	_mkdir(current_work_dir);*/
	//-------------------------------------------------
	string filename = "Result/coefficent.plt";
	//string filepress = "Result/eps_pressure.plt";
	string filelog = "Result/log.txt";
	log.open(filelog, ios::out);
	fDrug.open("Result/forceDrug.plt");
	fLift.open("Result/forceLift.plt");
	SetLog(log, grid, m, Re, alpha_f, beta_f, Zeidel_eps);
	log << endl;


	ApplyInitialData(U_new, grid); // Applying initial data to velocity 
	U_n = U_new;
	U_prev = U_new;



	cout << "flow over cylinder has been started\n";
	CalculateForce_X(Force_x, solidList, U_new, Cd, grid, alpha_f, beta_f, m);
	CalculateForce_Y(Force_y, solidList, V_new, Cl, grid, alpha_f, beta_f, m);
	fDrug << n << " " << Summ(Force_x) << endl;
	fLift << n << " " << Summ(Force_y) << endl;
	
	OutputVelocity_U(U_new, -1, output_step, solidList, grid);
	OutputVelocity_V(V_new, -1, output_step, solidList, grid);

	while (n <= N_max) {
		eps_u = 0.0;
		eps_v = 0.0;

		//<---------- prediction of velocity --------------------------
		B_u = CalculateB_u(U_n, V_n, U_prev, V_prev, P, Force_x, grid, Re);
		B_v = CalculateB_v(U_n, V_n, U_prev, V_prev, P, Force_y, grid, Re);
#pragma omp parallel sections num_threads(2)
		{

#pragma omp section
			{
				BiCGStab(U_new, grid.N1, grid.N2 + 1, OperatorA_u, B_u, grid,true);
			}
#pragma omp section
			{
				BiCGStab(V_new, grid.N1 + 1, grid.N2, OperatorA_v, B_v, grid,true);
			}

		}
		//ExplicPredVel(U_new,V_new,U_n,V_n,P,Force_x,Force_y,grid);

		//<----------end of prediction of velocity --------------------

		P_Right = Calculate_Press_Right(U_n, V_n, grid);
		for (int i = 0; i < (int)Delta_P.size(); ++i)
			for (int j = 0; j < (int)Delta_P[i].size(); ++j) {
				Delta_P[i][j] = 0.0;
			}

		

		eps_p = Calculate_Press_correction(Delta_P, P_Right, N_Zeidel, Zeidel_eps, grid,true);
		press_output << n << ' ' << eps_p << endl; //<---writing in closed stream



		for (int i = 0; i < grid.N1 + 1; ++i) {
			for (int j = 0; j < grid.N2 + 1; ++j) {
				P[i][j] = P[i][j] + 0.8 * Delta_P[i][j];
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

		//------------calculaating eps_u--------------------------
		for (int i = 0; i < grid.N1; ++i) {
			for (int j = 0; j < grid.N2 + 1; ++j) {
				if (fabs(U_n[i][j] - U_new[i][j]) > eps_u) {
					eps_u = fabs(U_n[i][j] - U_new[i][j]);
				}

				U_prev[i][j] = U_n[i][j];
				U_n[i][j] = U_new[i][j];
			}
		}
		//------------calculaating eps_v--------------------------
		for (int i = 0; i < grid.N1 + 1; ++i) {
			for (int j = 0; j < grid.N2; ++j) {
				if (fabs(V_n[i][j] - V_new[i][j]) > eps_v) {
					eps_v = fabs(V_n[i][j] - V_new[i][j]);
				}
				V_prev[i][j] = V_n[i][j];
				V_n[i][j] = V_new[i][j];
			}
		}
		//--------------------------------------------------------

		CalculateForce_X(Force_x, solidList, U_new, Cd, grid, alpha_f, beta_f, m);
		CalculateForce_Y(Force_y, solidList, V_new, Cl, grid, alpha_f, beta_f, m);
		
		fDrug << n + 1 << " " << Summ(Force_x) << endl;
		fLift << n + 1 << " " << Summ(Force_y) << endl;
		CreateMatrix(F_x, grid.N1, grid.N2 + 1);
		CreateMatrix(F_y, grid.N1 + 1, grid.N2);
		F_x = Calculate_F_real_x(U_n, V_n, U_prev, P, grid, Re);
		F_y = Calculate_F_real_y(U_n, V_n, V_prev, P, grid, Re);


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


	for (int i = 0; i < grid.N1; ++i) {
		for (int j = 0; j < grid.N2+1; ++j) {
			u[i][j] = 1.0;
		}
	}
}

int sgn(double x)
{
	(x >= 0) ? x = 1 : x = -1;
	return x;
}

Circle InputData(Grid& grid, double &M, int &Re, double &alpha_f, double &beta_f, double& Zeidel_eps, int& output_step, int& N_max, int& N_Zeidel) {

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
	double x, y, r;
	input >> x;
	input >> y;
	input >> r;
	Circle cylinder(x, y, r, 0, grid);
	input.close();


	grid.d_x = grid.L / (grid.N1 - 1);
	grid.d_y = grid.H / (grid.N2 - 1);

	return cylinder;
}

double Summ(Matrix& force) {
	double sum = 0;
	for (int i = 0; i < (int)force.size(); i++)
		for (int j = 0; j < (int)force[0].size(); j++)
		{
			sum += force[i][j];
		}


	return abs(sum);
}
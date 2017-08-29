#include "stdafx.h"
#include "SolidBody.h"
#include "Parameters.h"
#include "CalculateForce.h"
#include "Output.h"
#include "BiCGStab.h"
#include "Calculate_press.h"
#include "PredictVel.h"



#pragma warning(disable : 4996)//for using <chrono>
#pragma warning(disable : 4244)//for GetInfluenceArea


bool Debug = false;
bool InelasticCollision = false; //Perfectly inelastic collision --- абсолютно неупругие столкновения

using namespace std;


// fuctions

void InputData(Param& par);
void SetLog(ostream &log, Param par);
void PushLog(ostream &log, int n, double eps_u, double eps_v);
void ApplyInitialData(Matrix& u, Param par);
double ux_Poiseuille(double y, double H);
int sgn(double x);



int main() {

	const double epsilon = 1e-3;

	// declaring variables
	Param par;
	double eps_u = 0.0;
	double eps_v = 0.0;
	double eps_p = 0.0;



	int n = 0; // iteration counter
	InputData(par); // Get value of some variables
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
	for (int i = 0; i < 5; i++) {
		OperatorA_u[i].resize(par.N1);
		for (int j = 0; j < par.N1; j++) {
			OperatorA_u[i][j].resize(par.N2 + 1);
			fill(OperatorA_u[i][j].begin(), OperatorA_u[i][j].end(), 0);
		}
	}

	Matrix OperatorA_v[5];
	for (int i = 0; i < 5; i++) {
		OperatorA_v[i].resize(par.N1 + 1);
		for (int j = 0; j < par.N1 + 1; j++) {
			OperatorA_v[i][j].resize(par.N2);
			fill(OperatorA_v[i][j].begin(), OperatorA_v[i][j].end(), 0);
		}
	}
	Calculate_A_u(OperatorA_u, par, par.Re);
	Calculate_A_v(OperatorA_v, par, par.Re);


	// list of immersed solids
	list<Circle> solidList;
	ofstream output; // for Drag and Lift coefficents
	ofstream press_output; // press
	ofstream log;
	//-----------creating Result folder --------------
	//char current_work_dir[FILENAME_MAX];
	//_getcwd(current_work_dir, sizeof(current_work_dir));
	//strcat_s(current_work_dir, "\\Result");
	//_mkdir(current_work_dir);
	//-------------------------------------------------
	string filename = "Result/coefficent.plt";
	//string filepress = "Result/eps_pressure.plt";
	string filelog = "Result/log.txt";
	log.open(filelog, ios::out);
	SetLog(log, par);
	log << endl;


	ApplyInitialData(U_new, par); // Applying initial data to velocity 
	U_n = U_new;
	U_prev = U_new;


	//Firstly adding some circles
	GeomVec uc;
	std::fill(uc.begin(), uc.end(), 0.0);

	uc[1] = ux_Poiseuille(2.1, par.H);
	Circle c1(3.5, 2.1, par.R, par.NF, uc);

	uc[1] = ux_Poiseuille(4.9, par.H);
	Circle c2(3.5, 4.9, par.R, par.NF, uc);

	uc[1] = ux_Poiseuille(1.9, par.H);
	Circle c3(1.5, 1.9, par.R, par.NF, uc);

	uc[1] = ux_Poiseuille(5.1, par.H);
	Circle c4(1.5, 5.1, par.R, par.NF, uc);

	solidList.push_back(c1);
	solidList.push_back(c2);
	solidList.push_back(c3);
	solidList.push_back(c4);

	CalculateForce(Force_x, Force_y, solidList, U_new, V_new, par);

	OutputVelocity_U(U_new, -1, solidList, par);
	OutputVelocity_V(V_new, -1, solidList, par);

	while (n <= par.N_max) {

		//creation new solids
		if (n > 0 && fmod(n*par.d_t, 1.5) == 0.0) {
			int chance = 80;
			int rnd;
			rnd = rand() % 100 + 1;
			if (rnd <= chance) {
				double x = 1 + ((rand() % 100 + 1) / 100.0);
				double y = 1 + ((rand() % 200 + 1) / 100.0);
				uc[1] = ux_Poiseuille(y, par.H);
				Circle c(x, y, par.R, par.NF, uc);
				solidList.push_back(c);
			}
			rnd = rand() % 100 + 1;
			if (rnd <= chance) {
				double x = 1 + ((rand() % 100 + 1) / 100.0);
				double y = 4 + ((rand() % 200 + 1) / 100.0);
				uc[1] = ux_Poiseuille(y, par.H);
				Circle c(x, y, par.R, par.NF, uc);
				solidList.push_back(c);
			}
			rnd = rand() % 100 + 1;
			if (rnd <= chance) {
				double x = 3 + ((rand() % 100 + 1) / 100.0);
				double y = 1 + ((rand() % 200 + 1) / 100.0);
				uc[1] = ux_Poiseuille(y, par.H);
				Circle c(x, y, par.R, par.NF, uc);
				solidList.push_back(c);
			}

			rnd = rand() % 100 + 1;
			if (rnd <= chance) {
				double x = 3 + ((rand() % 100 + 1) / 100.0);
				double y = 4 + ((rand() % 200 + 1) / 100.0);
				uc[1] = ux_Poiseuille(y, par.H);
				Circle c(x, y, par.R, par.NF, uc);
				solidList.push_back(c);
			}
		}

		eps_u = 0.0;
		eps_v = 0.0;

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

		for (int i = 0; i < (int)Delta_P.size(); ++i)
			for (int j = 0; j < (int)Delta_P[i].size(); ++j) {
				Delta_P[i][j] = 0.0;
			}


		eps_p = Calculate_Press_correction(Delta_P, P_Right, par,false);

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

		//------------calculaating eps_u--------------------------
		for (int i = 0; i < par.N1; ++i) {
			for (int j = 0; j < par.N2 + 1; ++j) {
				if (fabs(U_n[i][j] - U_new[i][j]) > eps_u) {
					eps_u = fabs(U_n[i][j] - U_new[i][j]);
				}

				U_prev[i][j] = U_n[i][j];
				U_n[i][j] = U_new[i][j];
			}
		}
		//------------calculaating eps_v--------------------------
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

		CalculateForce(Force_x, Force_y, solidList, U_new, V_new, par);



		//--------------COLLISION CHECK---------------------------
		list<Circle>::iterator first;
		list<Circle>::iterator second;
		///--------------collisions between particles---------------------------------
		for (auto one = solidList.begin(); one != solidList.end(); one++) {
			for (auto two = next(one); two != solidList.end(); two++) {
				double distance = length(one->xc - two->xc); //<----distance between two particles 
				if (one->xc[1] < two->xc[1]) {
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
						first ->uc = (first->uc + second->uc) / 2;
						second->uc =  first->uc;
					}
					else {
						//this is perfectly elastic impact
						double v1, v2;
						v1 = sgn(first ->uc[1])*length(first ->uc);
						v2 = sgn(second->uc[1])*length(second->uc);

						double tetta_1, tetta_2, phi;
						if (v1 != 0) { tetta_1 = atan(first->uc[2] / first->uc[1]); }
						else { tetta_1 = 0; } //??if (v1 == 0) { tetta_1 = M_PI_2 - tetta_2 / 2; }
						if (v2 != 0) { tetta_2 = atan(second->uc[2] / second->uc[1]); }
						else { tetta_2 = 0; } //??if (v2 == 0) { tetta_2 = M_PI_2 - tetta_1 / 2; }
						phi = atan((second->xc[2] - first->xc[2]) / (second->xc[1] - first->xc[1]));
						first->uc[1] = v2*cos(tetta_2 - phi)*cos(phi) + v1*sin(tetta_1 - phi)*cos(phi + M_PI_2);
						first->uc[2] = v2*cos(tetta_2 - phi)*sin(phi) + v1*sin(tetta_1 - phi)*sin(phi + M_PI_2);
						second->uc[1] = v1*cos(tetta_1 - phi)*cos(phi) + v2*sin(tetta_2 - phi)*cos(phi + M_PI_2);
						second->uc[2] = v1*cos(tetta_1 - phi)*sin(phi) + v2*sin(tetta_2 - phi)*sin(phi + M_PI_2);
					}
				}
			}
		}
		///--------------end of collisions between particles---------------------------------

		///-------------collision with walls--------------------------
		for (auto one = solidList.begin(); one != solidList.end(); one++) {
			double DistUpper = par.H - one->xc[2];//<----distance to upper wall
			double DistLower = one->xc[2];//<-------distance to lower wall
			if (DistUpper < one->r || DistLower < one->r) {
				if (Debug) cout << "COLLISION DETECTED";
				one->uc[2] = -one->uc[2];
			}
		}
		///-------------end of collision with walls--------------------------

//--------------END OF COLLISION CHECK---------------------------------------------------


		for (auto it = solidList.begin(); it != solidList.end();) {
			if ((it->moveSolid == false)) {
				it->moveSolid = true;
			}

			if (it->moveSolid) {
				//update position
				for (int k = 0; k < par.NF; ++k) {
					//rotate
					GeomVec r = it->Nodes[k].x - it->xc;
					GeomVec x_temp = rotate_Vector_around_vector(r, it->omega  * length(r) * par.d_t); //
					it->Nodes[k].x = it->xc + x_temp; // rotate solid by angle $omega$ * $dt$
					// move
					it->Nodes[k].x += it->uc * par.d_t;
				}
				it->xc += it->uc * par.d_t;
			}

			//delete bodies which move 95% of length
			if (it->xc[1] > par.L*0.95) {
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

		if (0 == n % par.output_step) {
			OutputVelocity_U(U_new, n, solidList, par);
			OutputVelocity_V(V_new, n, solidList, par);
			OutputPressure(P, n, solidList, par);
		}



		if (eps_u < epsilon && eps_v < epsilon) {
			OutputVelocity_U(U_new, n, solidList, par);
			OutputVelocity_V(V_new, n, solidList, par);
			OutputPressure(P, n, solidList, par);
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
void SetLog(ostream& log, Param par) {

	log << "The IBM program starts.		";
	time_t t = chrono::system_clock::to_time_t(chrono::system_clock::now());   // get time now
	log << ctime(&t) << endl;
	log << "The parameters are the following:" << endl;
	log << "Reynolds number               : Re  = " << par.Re << endl;
	log << "Channel length                : L   = " << par.L << endl;
	log << "Channel width                 : W   = " << par.H << endl;
	log << "Number of nodes on            : N1  = " << par.N1 << endl;
	log << "Number of nodes on            : N2  = " << par.N2 << endl;
	log << "Number of nodes for a particle: NF  = " << par.NF << endl;
	log << "Time step                     : tau = " << par.d_t << endl;
	log << "Force parameter alpha         : alpha = " << par.alpha_f << endl;
	log << "Force parameter beta          : beta  = " << par.beta_f << endl;
	log << "Tolerance for Zeidel method   : tol = " << par.Zeidel_eps << endl;

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
void ApplyInitialData(Matrix &u, Param par) {

	// Poiseuille flow 
	for (int i = 0; i < par.N1; ++i) {
		for (int j = 1; j < par.N2; ++j) {
			u[i][j] = ux_Poiseuille((j - 0.5)*par.d_y, par.H);
		}
	}
}

double ux_Poiseuille(double y, double H) {
	double ux = (pow(H / 2.0, 2) - pow(y - H / 2.0, 2));
	return ux;
}

int sgn(double x)
{
	(x >= 0) ? x = 1 : x = -1;
	return x;
}


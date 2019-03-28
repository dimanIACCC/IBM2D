#include "stdafx.h"
#include "BodyOfProgram.h"

void BodyOfProgram(Param par, std::list<Circle> solidList, Matrix U_n, Matrix V_n, Matrix P_n, bool NeedNewLog) {
																		 
#pragma region SetMatrices 
	CreateMatrix(U_new, par.N1_u, par.N2_u);
	CreateMatrix(V_new, par.N1_v, par.N2_v);
	CreateMatrix(P_new, par.N1_p, par.N2_p);
	CreateMatrix(Fx_n  , par.N1_u, par.N2_u);
	CreateMatrix(Fy_n  , par.N1_v, par.N2_v);
	CreateMatrix(Fx_new, par.N1_u, par.N2_u);
	CreateMatrix(Fy_new, par.N1_v, par.N2_v);

	CreateMatrix(U_exact, par.N1_u, par.N2_u);
	CreateMatrix(V_exact, par.N1_v, par.N2_v);
	CreateMatrix(P_exact, par.N1_p, par.N2_p);
	CreateMatrix(U_exact_n, par.N1_u, par.N2_u);
	CreateMatrix(V_exact_n, par.N1_v, par.N2_v);
	CreateMatrix(P_exact_n, par.N1_p, par.N2_p);
	CreateMatrix(U      , par.N1_u, par.N2_u);
	CreateMatrix(V      , par.N1_v, par.N2_v);
	CreateMatrix(P      , par.N1_p, par.N2_p);
	CreateMatrix(P_d    , par.N1_p, par.N2_p);
	CreateMatrix(P_old  , par.N1_p, par.N2_p);
	CreateMatrix(dU, par.N1_u, par.N2_u);
	CreateMatrix(dV, par.N1_v, par.N2_v);
	CreateMatrix(dP, par.N1_p, par.N2_p);

	Template A_u;
	Template A_v;
	Calculate_A(A_u, par, par.Re);
	Calculate_A(A_v, par, par.Re);
#pragma endregion SetMatrices


	std::ofstream log;														// creation output stream for writing log.
	if (NeedNewLog) {														// condition for creation new log file or
		log.open(par.WorkDir + "log.txt", std::ios::out);					// continue writing if it`s after hibernation 
		SetLog(log, par);													// 
	}else
		log.open(par.WorkDir + "log.txt", std::ios::app);

	std::ofstream history;
	std::string filename = par.WorkDir + "/history.plt";
	if (par.BC == Taylor_Green || par.BC == Lamb_Oseen || par.BC == Line_Vortex) {
		history.open(filename);
		history << "title = history" << std::endl;
		history << "Variables = n error_P error_U error_V" << std::endl;
	}

	for ( ; par.N_step <= par.N_max; ++par.N_step) {                                 // main cycle of time iterations
		if (par.N_step % 1000 == 0)
			MakeHibernationFile(par, solidList, U_n, V_n, P_n);              // writting hibernation file for prior time step
		if (par.N_step % par.output_step == 0 || par.N_step < 1000)
			Output(P_n, U_n, V_n, Fx_n, Fy_n, par.N_step, solidList, par);

		Add_Solids(solidList, par);                                                     // add solids if the conditions are fulfilled

		//if (par.N_step == 1000) {
		//	par.d_t /= 10.0; //workaround
		//	Calculate_A(A_u, par, par.Re);
		//	Calculate_A(A_v, par, par.Re);
		//}

		Calculate_u_p(U_n, U_new, V_n, V_new, P_n, P_new, Fx_n, Fx_new, Fy_n, Fy_new, A_u, A_v, solidList, par);  // calculate velocity and pressure at the new time step

		// Output_eq_terms("eq_terms", par.N_step, U_n, V_n, U_new, V_new, P_n, P_new, Fx_new, par, Du);

		double eps_u = diff(U_n, U_new);												// calculate maximum difference between current and prior velocity fields
		double eps_v = diff(V_n, V_new);

		// single averaging
		U = (U_n + U_new) * 0.5;
		V = (V_n + V_new) * 0.5;
		P = (P_n + P_new) * 0.5;

		// step (n + 0.5)
		fill_exact(U_exact, V_exact, P_exact, par, par.d_t*(par.N_step + 0.5));
		dU = U - U_exact;
		dV = V - V_exact;
		dP = P - P_exact;
		double max_dU = max(dU);
		double max_dV = max(dV);
		double max_dP = max(dP);

		// double averaging
		if (par.N_step == 0) P_d = P_n;
		else                 P_d = (P_old + P) * 0.5;
		P_old = P;

		//step n
		fill_exact(U_exact_n, V_exact_n, P_exact_n, par, par.d_t*(par.N_step));
		double max_dU_n = max(U_n - U_exact_n);
		double max_dV_n = max(V_n - V_exact_n);
		double max_dP_d = max(P_d - P_exact_n);

		history << par.N_step << "  " << max_dP / max(P_exact) << "  " << max_dU_n / max(U_exact_n) << "  " << max_dV_n / max(V_exact_n) << std::endl;

		if (par.N_step % par.output_step == 0 || par.N_step < 1000) {
			//Output_U(dU, "dU", par.N_step, par);
			//Output_P(dP, "dP", par.N_step, par);

			//Output_P(P      , "P"      , par.N_step, par);
			//Output_P(P_exact, "P_exact", par.N_step, par);
			//Output_P(P_n    , "P_n"    , par.N_step, par);
			//Output_P(P_new  , "P_new"  , par.N_step, par);

			//Output_U(U      , "U"      , par.N_step, par);
			//Output_U(U_exact, "U_exact", par.N_step, par);
			//Output_U(U_n    , "U_n"    , par.N_step, par);
			//Output_U(U_new  , "U_new"  , par.N_step, par);

			//Output(P, U, V, Fx_new, Fy_new, par.N_step, solidList, par);
			//Output_U(U, "U", par.N_step, par);
			//Output_V(V, "V", par.N_step, par);

			//Output_U(Fx_new, "U", par.N_step, par);
			//Output_V(Fy_new, "V", par.N_step, par);
		}

		U_n = U_new;
		V_n = V_new;
		P_n = P_new;
		
		//workaround for moving walls
		//for (auto& solid : solidList) {
		//	for (size_t i = 0; i < U_n.size(); i++)
		//		for (size_t j = 0; j < U_n[0].size(); j++)
		//			U_n[i][j] = U_n[i][j] - solid.u[1];
		//			par.u_wall -= solid.u[1];
		//			solid.u[1] = 0;
		//}

		Solids_move(solidList, par);												// moving solids if it is necessary (checking it up inside)
																					// and detection of collisions

		PushLog(log, par.N_step, eps_u, eps_v);                                              // writting log into log-file
		log.flush();


		const double epsilon = 1e-7;
		if (eps_u < epsilon && eps_v < epsilon && par.N_step > 1000) {
			break;
		}
	}

	
	log.close();
}

void MakeHibernationFile(Param& par, std::list<Circle>& solidList, Matrix& U_n, Matrix& V_n, Matrix& P_n) {
	std::ofstream output;
	std::string filename = par.WorkDir + "step" + std::to_string(par.N_step) + ".txt";
	output.open(filename);
	output << std::setprecision(15);
	output << "par{" << std::endl;
	output << "N_step = " << par.N_step << std::endl;
	output << "Re = " << par.Re << std::endl;
	output << "L = " << par.L << std::endl;
	output << "H = " << par.H << std::endl;
	output << "N1 = " << par.N1 << std::endl;
	output << "N2 = " << par.N2 << std::endl;
	output << "d_t = " << par.d_t << std::endl;
	output << "Nn = " << par.Nn << std::endl;
	output << "rho = " << par.rho << std::endl;
	output << "r = " << par.r << std::endl;
	output << "output_step = " << par.output_step << std::endl;
	output << "N_max = " << par.N_max << std::endl;
	output << "DeltaP_method = " << par.DeltaP_method << std::endl;
	output << "N_Zeidel = " << par.N_Zeidel << std::endl;
	output << "Zeidel_eps = " << par.Zeidel_eps << std::endl;
	output << "eps_P = " << par.eps_P << std::endl;
	output << "InelasticCollision = " << par.InelasticCollision << std::endl;
	output << "k_dist = " << par.k_dist << std::endl;
	output << "AddSolids_N = " << par.AddSolids_N << std::endl;
	output << "AddSolids_start = " << par.AddSolids_start << std::endl;
	output << "AddSolids_interval = " << par.AddSolids_interval << std::endl;
	output << "BC = " << par.BC << std::endl;
	output << "output_step = " << par.output_step << std::endl;
	output << "N_Force = " << par.N_Force << std::endl;
	output << "SolidName_max = " << par.SolidName_max << std::endl;
	output << "}" << std::endl;

	output << "U_n{" << std::endl;
	for (int i = 0; i < par.N1_u; i++) {
		for (int j = 0; j < par.N2_u; j++) output << U_n[i][j] << ' ';
		output << std::endl;
	}
	output << "}" << std::endl;

	output << "V_n{" << std::endl;
	for (int i = 0; i < par.N1_v; i++) {
		for (int j = 0; j < par.N2_v; j++) output << V_n[i][j] << ' ';
		output << std::endl;
	}
	output << "}" << std::endl;

	output << "P_n{" << std::endl;
	for (int i = 0; i < par.N1_p; i++) {
		for (int j = 0; j < par.N2_p; j++) output << P_n[i][j] << ' ';
		output << std::endl;
	}
	output << "}" << std::endl;

	output << "<Solidlist>" << std::endl;
	for (auto one = solidList.begin(); one != solidList.end(); one++) {
		output << "<Solid>" << std::endl;
		output << "moving = " << one->moving << std::endl;
		output << "x_n = " << std::endl << one->x_n << std::endl;
		output << "u_n = " << std::endl << one->u_n << std::endl;
		output << "omega_n = " << std::endl << one->omega_n << std::endl;
		output << "I = " << one->I << std::endl;
		output << "rho = " << one->rho << std::endl;
		output << "V = " << one->V << std::endl;
		output << "Nn = " << one->Nn << std::endl;
		output << "name = " << one->name << std::endl;
		output << "r = " << one->r << std::endl;
		output << "n_moving" << one->n_moving << std::endl;
		output << "<Nodes>" << std::endl;
		for (int j = 0; j < par.Nn; j++) {
			output << "Node{" << std::endl;
			output << "x_n = " << std::endl << one->Nodes[j].x_n << std::endl;
			output << "n = " << std::endl << one->Nodes[j].n << std::endl;
			output << "}" << std::endl;
		}
		output << "<\\Nodes>" << std::endl;
		output << "<\\Solid>" << std::endl;
	}
	output << "<\\Solidlist>" << std::endl;
	output.close();
}


//this method is restoring calculated values to continue calculations
void Awake(std::string &filename, Param &par, std::list<Circle>& solidList, Matrix& U_n, Matrix& V_n, Matrix& P_n) {
	std::ifstream hibernation_source;
	std::string line;
	std::string PAR, VALUE;
	char ch;


	std::cout << filename << std::endl;
	hibernation_source.open(filename);
	if (hibernation_source.is_open()) {
		while (getline(hibernation_source, line)) { // read line from file to string $line$
			GetParValue(line, PAR, VALUE);
			if (line == "par{") {
				while (line != "}") {
					getline(hibernation_source, line);
					if (line == "}") break;
					par.read_line(line);
				}
				par.init();
			}
			else if (line == "U_n{") {
				U_n.resize(par.N1_u);
				for (int i = 0; i < par.N1_u; i++) U_n[i].resize(par.N2_u);
				for (int i = 0; i < par.N1_u; i++)
					for (int j = 0; j < par.N2_u; j++) hibernation_source >> U_n[i][j];


				if (hibernation_source >> ch && ch != '}') std::cout << "Some troubles in U_n" << std::endl;
			}
			else if (line == "V_n{") {
				V_n.resize(par.N1_v);
				for (int i = 0; i < par.N1_v; i++) V_n[i].resize(par.N2_v);
				for (int i = 0; i < par.N1_v; i++)
					for (int j = 0; j < par.N2_v; j++) hibernation_source >> V_n[i][j];

				if (hibernation_source >> ch && ch != '}') std::cout << "Some troubles in V_n" << std::endl;
			}
			else if (line == "P_n{") {
				P_n.resize(par.N1_p);
				for (int i = 0; i < par.N1_p; i++) P_n[i].resize(par.N2_p);
				for (int i = 0; i < par.N1_p; i++)
					for (int j = 0; j < par.N2_p; j++) hibernation_source >> P_n[i][j];

				if (hibernation_source >> ch && ch != '}') std::cout << "Some troubles in P" << std::endl;
			}
			else if (line == "<Solidlist>")
				while (line != "<\\Solidlist>") {
					getline(hibernation_source, line);
					Circle c(0, 0, par);
					if (c.name > par.SolidName_max) par.SolidName_max = c.name;
					if (line == "<\\Solidlist>") break;
					if (line == "<Solid>") {
						while (line != "<\\Solid>") {
							getline(hibernation_source, line);
							if (line == "<\\Solid>") {
								solidList.push_back(c);
								break;
							}
							double x = par.L*0.1;
							double y = par.H*0.5;
							double ux = 0;
							double uy = 0;
							double omega = 0;
							double rho = par.rho;
							int Nn = par.Nn;
							bool moving = true;
							double r = par.r;

							GetParValue(line, PAR, VALUE);

							if      (PAR == "moving")        c.moving = bool(stoi(VALUE));
							else if (PAR == "x_n")           hibernation_source >> c.x_n;//
							else if (PAR == "u_n")           hibernation_source >> c.u_n;//
							else if (PAR == "omega_n")       hibernation_source >> c.omega_n;//
							else if (PAR == "I")             c.I = stod(VALUE);
							else if (PAR == "rho")           c.rho = stod(VALUE);
							else if (PAR == "V")             c.V = stod(VALUE);
							else if (PAR == "Nn")            c.Nn = stoi(VALUE);
							else if (PAR == "name")          c.name = stoi(VALUE);
							else if (PAR == "r")             c.r = stod(VALUE);
							else if (PAR == "n_moving")      c.n_moving = stoi(VALUE);
							else if (PAR == "<Nodes>") {
								while (line != "<\\Nodes>") {
									getline(hibernation_source, line);
									if (line == "<\\Nodes>") break;
									for (int j = 0; j < par.Nn; j++) {
										if (line == "Node{") {
											while (line != "}") {
												getline(hibernation_source, line);
												if (line == "}") {
													getline(hibernation_source, line);
													break;
												}
												GetParValue(line, PAR, VALUE);
												if (PAR == "x_n")				hibernation_source >> c.Nodes[j].x_n;
												else if (PAR == "n")			hibernation_source >> c.Nodes[j].n;
												c.Nodes[j].x = c.Nodes[j].x_n;
											}
										}
									}
								}
							}
							c.x = c.x_n;
							c.u = c.u_n;
							c.omega = c.omega_n;
						}
					}
				}
		}
	}
	else std::cout << "Hibernation file is not found. Please check the path\n The correct command is like -h=ResultFolder" << std::endl;
	std::cout << "End of awaking" << std::endl;

}
#include "stdafx.h"
#include "BodyOfProgram.h"

void BodyOfProgram(Param &par, std::vector<Solid> &solidList, std::vector<Node> &Nodes, Matrix &U_n, Matrix &V_n, Matrix &P) {
																		 
#pragma region SetMatrices 
	CreateMatrix(U_new, par.N1_u, par.N2_u);
	CreateMatrix(V_new, par.N1_v, par.N2_v);
	CreateMatrix(Fx   , par.N1_u, par.N2_u);
	CreateMatrix(Fy   , par.N1_v, par.N2_v);

	CreateMatrix(U_exact, par.N1_u, par.N2_u);
	CreateMatrix(V_exact, par.N1_v, par.N2_v);
	CreateMatrix(P_exact, par.N1_p, par.N2_p);

	CreateMatrix(dU, par.N1_u, par.N2_u);
	CreateMatrix(dV, par.N1_v, par.N2_v);
	CreateMatrix(dP, par.N1_p, par.N2_p);

#pragma endregion SetMatrices

	std::ofstream history;
	std::string filename = par.WorkDir + "/history.plt";
	history_init(par.WorkDir, "history", par.BC);

	for ( ; par.N_step <= 5000000; ++par.N_step) {                                 // main cycle of time iterations

		Add_Solids(solidList, Nodes, par);                                                     // add solids if the conditions are fulfilled
		std::cout << "Add_Solids finished" << std::endl;
		Output(P, U_new, V_new, Fx, Fy, par.N_step, solidList, Nodes, par);
		std::getchar();

		if (par.N_step % 100 == 0)
			MakeHibernationFile(par, solidList, Nodes, U_n, V_n, P);              // writting hibernation file for prior time step

		Calculate_u_p(U_n, U_new, V_n, V_new, P, Fx, Fy, solidList, Nodes, par);  // calculate velocity and pressure at the new time step

		Solids_position_new(solidList, Nodes, par);

		// Output_eq_terms("eq_terms", par.N_step, U_n, V_n, U_new, V_new, P_n, P_new, Fx_new, par, Du);

		double eps_u = diff(U_n, U_new);												// calculate maximum difference between current and prior velocity fields
		double eps_v = diff(V_n, V_new);

		if (par.N_step % par.output_step == 0 || par.N_step < 1)
			Output(P, U_new, V_new, Fx, Fy, par.N_step, solidList, Nodes, par);

		fill_exact(U_exact, V_exact, P_exact, par, par.d_t*(par.N_step + 1), par.d_t*(par.N_step + 0.5));
		dU = U_new - U_exact;
		dV = V_new - V_exact;
		dP = P     - P_exact;
		double max_dU = Matrix_max(dU);
		double max_dV = Matrix_max(dV);
		double max_dP = Matrix_max(dP);

		if (par.N_step % par.output_step == 0 || par.N_step < 1000) {
			//Output_U(dU, "dU", par.N_step, par);
			//Output_P(dP, "dP", par.N_step, par);

			//Output(P, U_new, V_new, Fx, Fy, par.N_step, solidList, par);
		}

		U_n = U_new;
		V_n = V_new;
		
		//workaround for moving walls
		//for (auto& solid : solidList) {
		//	for (size_t i = 0; i < U_n.size(); i++)
		//		for (size_t j = 0; j < U_n[0].size(); j++)
		//			U_n[i][j] = U_n[i][j] - solid.u[1];
		//			par.u_wall -= solid.u[1];
		//			solid.u[1] = 0;
		//}


		Solids_move(solidList, Nodes, par);												// moving solids if it is necessary (checking it up inside)
																					        // and detection of collisions
		double h_average=0;
		//h_average_of_Solids_Layer(solidList, par, h_average);
		
		if (par.BC == box) {
			history_log(par.WorkDir, "history", par.N_step*par.d_t, h_average, 0., 0.);
		}
		else if (par.BC == Taylor_Green || par.BC == Lamb_Oseen || par.BC == Line_Vortex) {
			history_log(par.WorkDir, "history", par.N_step*par.d_t, max_dP / Matrix_max(P_exact), max_dU / Matrix_max(U_exact), max_dV / Matrix_max(V_exact) );
		}

		std::cout << "n = " << std::setw(6) << par.N_step << std::endl;

		par.time += par.d_t;

		const double epsilon = 1e-13;
		if (eps_u < epsilon && eps_v < epsilon && par.N_step > 1000) {
			std::cout << "u=0" << std::endl;
			break;
		}
	}

	
	//log.close();
}

void MakeHibernationFile(Param& par, std::vector<Solid>& solidList, std::vector<Node> &Nodes, Matrix& U_n, Matrix& V_n, Matrix& P_n) {
	std::ofstream output;
	std::string filename = par.WorkDir + "step" + std::to_string(par.N_step) + ".txt";
	output.open(filename);
	output << std::setprecision(15);
	output << "par{" << std::endl;

	output << "N_step = " << par.N_step << std::endl;
	output << "time = " << par.time << std::endl;
	output << "BC = " << par.BC << std::endl;

	// physical parameters
	output << "d_t = " << par.d_t << std::endl;
	output << "L = " << par.L << std::endl;
	output << "H = " << par.H << std::endl;
	output << "Re = " << par.Re << std::endl;
	output << "grad_p_x = " << par.grad_p_x << std::endl;
	output << "Gravity_angle = " << par.Gravity_angle *180. / M_PI << std::endl;
	output << "Gravity_module = " << par.Gravity_module << std::endl;

	// numerical parameters
	output << "N1 = " << par.N1 << std::endl;
	output << "N2 = " << par.N2 << std::endl;
	output << "Nn = " << par.Nn << std::endl;
	output << "output_step = " << par.output_step << std::endl;
	output << "IBM = " << par.IBM << std::endl;
	output << "DeltaP_method = " << par.DeltaP_method << std::endl;
	output << "N_Zeidel = " << par.N_Zeidel << std::endl;
	output << "Zeidel_eps = " << par.Zeidel_eps << std::endl;
	output << "s_max = " << par.s_max << std::endl;
	output << "eps_P = " << par.eps_P << std::endl;

	// parameters for many particles
	output << "rho = " << par.rho << std::endl;
	output << "r = " << par.r << std::endl;
	output << "AddSolids_N = " << par.AddSolids_N << std::endl;
	output << "AddSolids_start = " << par.AddSolids_start << std::endl;
	output << "AddSolids_interval = " << par.AddSolids_interval << std::endl;
	output << "SolidName_max = " << par.SolidName_max << std::endl;
	output << "k_dist = " << par.k_dist << std::endl;

	// parameters for special problems
	output << "Lamb_Oseen_r0 = " << par.Lamb_Oseen_r0 << std::endl;
	output << "u_in        = " << par.u_in << std::endl;
	output << "u_down = " << par.u_down << std::endl;
	output << "u_up   = " << par.u_up   << std::endl;
	output << "omega_BC = " << par.omega_BC << std::endl;

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
		output << "alpha = " << std::endl << one->alpha << std::endl;
		output << "I = " << one->I << std::endl;
		output << "rho = " << one->rho << std::endl;
		output << "V = " << one->V << std::endl;
		output << "Nn = " << one->Nn << std::endl;
		output << "name = " << one->name << std::endl;
		output << "r = " << one->r << std::endl;
		output << "<Nodes>" << std::endl;
		for (int j = 0; j < one->Nn; j++) {
			output << "Node{" << std::endl;
			output << "x_n = " << std::endl << Nodes[one->IndNodes[j]].x_n << std::endl;
			output << "ds = "  << std::endl << Nodes[one->IndNodes[j]].ds << std::endl;
			output << "}" << std::endl;
		}
		output << "<\\Nodes>" << std::endl;
		output << "<\\Solid>" << std::endl;
	}
	output << "<\\Solidlist>" << std::endl;
	output.close();
}


//this method is restoring calculated values to continue calculations
void Awake(std::string &filename, Param &par, std::vector<Solid>& solidList, std::vector<Node> &Nodes, Matrix& U_n, Matrix& V_n, Matrix& P_n) {
	std::ifstream hibernation_source;
	std::string line;
	std::string PAR, VALUE;
	char ch;

	std::cout << "Awake" << std::endl;

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

					if (line == "<\\Solidlist>") break;
					if (line == "<Solid>") {
						double x = par.L*0.1;
						double y = par.H*0.5;
						double ux = 0;
						double uy = 0;
						double omega = 0;

						GeomVec x_n, u_n, omega_n, alpha;

						double rho = par.rho;
						int name = par.SolidName_max+1;
						int Nn = par.Nn;
						int moving = 1;
						double shape = par.shape;
						double r = par.r;
						double e = par.e;
						bool Poiseuille = false;   //key for initial ux, uy and omega_new corresponding to Poiseuille flow

						while (line != "<\\Solid>") {
							getline(hibernation_source, line);
							GetParValue(line, PAR, VALUE);

							if      (PAR == "moving")        moving = stoi(VALUE);
							else if (PAR == "x_n")           hibernation_source >> x_n;//
							else if (PAR == "u_n")           hibernation_source >> u_n;//
							else if (PAR == "omega_n")       hibernation_source >> omega_n;//
							else if (PAR == "alpha")         hibernation_source >> alpha;
							else if (PAR == "rho")           rho = stod(VALUE);
							else if (PAR == "Nn")            Nn = stoi(VALUE);
							else if (PAR == "name")          name = stoi(VALUE);
							else if (PAR == "shape")         shape = stoi(VALUE);
							else if (PAR == "r")             r = stod(VALUE);
							else if (PAR == "e")             e = stod(VALUE);
							else if (PAR == "<Nodes>") {
								Nodes.resize(par.Nn_max + Nn);
								while (line != "<\\Nodes>") {
									getline(hibernation_source, line);
									if (line == "<\\Nodes>") break;
									for (int j = 0; j < Nn; j++) {
										int Ind = par.Nn_max + j;
										if (line == "Node{") {
											while (line != "}") {
												getline(hibernation_source, line);
												if (line == "}") {
													getline(hibernation_source, line);
													break;
												}
												GetParValue(line, PAR, VALUE);
												if (PAR == "x_n")               hibernation_source >> Nodes[Ind].x_n;
												if (PAR == "ds")                hibernation_source >> Nodes[Ind].ds;
												Nodes[Ind].x = Nodes[Ind].x_n;
											}
										}
									}
								}
							}
							if (line == "<\\Solid>") {
								double alpha0 = alpha[3];
								Solid c(x, y, ux, uy, alpha0, omega, rho, Nn, moving, name, shape, r, e);
								c.add_Nodes(Nodes, par.Nn_max);
								fill_solid_ds(Nodes, par.Nn_max, c.Nn, c.shape, 0.5*(par.d_x + par.d_y));
								par.Nn_max += c.Nn;

								if (c.name > par.SolidName_max) par.SolidName_max = c.name;

								c.x_n = x_n;
								c.x   = x_n;
								c.x_n_plt = x_n;
								c.x_plt   = x_n;
								c.u_n = u_n;
								c.u   = u_n;
								c.omega_n = omega_n;
								c.omega = omega_n;

								solidList.push_back(c);
								break;
							}
						}
					}
				}
		}
	}
	else std::cout << "Hibernation file is not found. Please check the path\n The correct command is like -h=ResultFolder" << std::endl;
	std::cout << "End of awaking" << std::endl;

}

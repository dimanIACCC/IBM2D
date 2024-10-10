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
	std::string filename = par.WorkDir + "/" + "history.plt";
	history_init(par.WorkDir, "history", par.BC);

	for ( ; par.N_step <= 5000000; ++par.N_step) {                                 // main cycle of time iterations

		Add_Solids(solidList, Nodes, par);                                                     // add solids if the conditions are fulfilled
		//Output(P, U_new, V_new, Fx, Fy, par.N_step, solidList, Nodes, par);
		//std::getchar();

		if (par.N_step % 10 == 0)
			Save_Data(par, solidList, Nodes, U_n, V_n, P);              // writting hibernation file for prior time step

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

}

void Save_Data(Param& par, std::vector<Solid>& solidList, std::vector<Node> &Nodes, Matrix& U_n, Matrix& V_n, Matrix& P_n) {
	std::ofstream output;
	std::string filename = par.WorkDir + "/" + "step" + std::to_string(par.N_step) + ".txt";
	output.open(filename);
	output << std::setprecision(15);
	output << "<par>" << std::endl;

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
	output << "Nn_ = " << par.Nn_ << std::endl;
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
	output << "r0 = " << par.r0 << std::endl;
	output << "e = " << par.e << std::endl;
	output << "AddSolids_N = " << par.AddSolids_N << std::endl;
	output << "AddSolids_start = " << par.AddSolids_start << std::endl;
	output << "AddSolids_interval = " << par.AddSolids_interval << std::endl;
	output << "SolidName_max = " << par.SolidName_max << std::endl;
	output << "k_u_dist = " << par.k_u_dist << std::endl;
	output << "k_r_dist = " << par.k_r_dist << std::endl;
	output << "k_u_collide = " << par.k_u_collide << std::endl;
	output << "k_r_collide = " << par.k_r_collide << std::endl;

	// parameters for special problems
	output << "Lamb_Oseen_r0 = " << par.Lamb_Oseen_r0 << std::endl;
	output << "u_in        = " << par.u_in << std::endl;
	output << "u_down = " << par.u_down << std::endl;
	output << "u_up   = " << par.u_up   << std::endl;
	output << "omega_BC = " << par.omega_BC << std::endl;

	output << "</par>" << std::endl;

	output << "<U>" << std::endl;
	for (int i = 0; i < par.N1_u; i++) {
		for (int j = 0; j < par.N2_u; j++) output << U_n[i][j] << ' ';
		output << std::endl;
	}
	output << "</U>" << std::endl;

	output << "<V>" << std::endl;
	for (int i = 0; i < par.N1_v; i++) {
		for (int j = 0; j < par.N2_v; j++) output << V_n[i][j] << ' ';
		output << std::endl;
	}
	output << "</V>" << std::endl;

	output << "<P>" << std::endl;
	for (int i = 0; i < par.N1_p; i++) {
		for (int j = 0; j < par.N2_p; j++) output << P_n[i][j] << ' ';
		output << std::endl;
	}
	output << "</P>" << std::endl;

	output << "<Solids>" << std::endl;
	for (auto one = solidList.begin(); one != solidList.end(); one++) {
		output << "<Solid>" << std::endl;
		output << "moving = " << one->moving << std::endl;
		output << "x = " << one->x_n[1] << std::endl;
		output << "y = " << one->x_n[2] << std::endl;
		output << "ux = " << one->u_n[1] << std::endl;
		output << "uy = " << one->u_n[2] << std::endl;
		output << "omega = " << one->omega_n[3] / M_PI * 180 << std::endl;
		output << "alpha = " << one->alpha[3]   / M_PI * 180 << std::endl;
		output << "rho = " << one->rho << std::endl;
		output << "Nn_r0 = " << one->Nn_r0 << std::endl;
		output << "Nn_r = " << one->Nn_r << std::endl;
		output << "Nn = " << one->Nn << std::endl;
		output << "name = " << one->name << std::endl;
		output << "shape = " << one->shape << std::endl;
		output << "r0 = " << one->r0 << std::endl;
		output << "r = " << one->r << std::endl;
		output << "e = " << one->e << std::endl;
		output << "<Nodes>" << std::endl;
		for (int j = 0; j < one->Nn; j++) {
			output << "<Node>" << std::endl;
			output << "x = "  << Nodes[one->IndNodes[j]].x_n[1] << std::endl;
			output << "y = "  << Nodes[one->IndNodes[j]].x_n[2] << std::endl;
			output << "ds = " << Nodes[one->IndNodes[j]].ds << std::endl;
			output << "</Node>" << std::endl;
		}
		output << "</Nodes>" << std::endl;
		output << "</Solid>" << std::endl;
	}
	output << "</Solids>" << std::endl;
	output.close();
}


//this method is restoring calculated values to continue calculations
void Load_Data(std::string &filename, Param &par, std::vector<Solid>& Solids, std::vector<Node> &Nodes, Matrix& U_n, Matrix& V_n, Matrix& P_n) {
	std::ifstream input;
	std::string line;
	bool key_U = FALSE;
	bool key_V = FALSE;
	bool key_P = FALSE;

	std::cout << "Loading Data from file = " << std::endl;
	std::cout << filename << std::endl;
	input.open(filename);
	if (input.is_open()) {
		while (getline(input, line)) { // read line from file to string $line$
			if (line == "<par>") {
				par.read(input);
			}
			else if (line == "<U>") {
				key_U = TRUE;
				Matrix_resize(U_n, par.N1_u, par.N2_u);
				Matrix_read  (U_n, par.N1_u, par.N2_u, input);
				if (input >> line && line != "</U>") std::cout << "Some troubles in <U>" << std::endl;
			}
			else if (line == "<V>") {
				key_V = TRUE;
				Matrix_resize(V_n, par.N1_v, par.N2_v);
				Matrix_read  (V_n, par.N1_v, par.N2_v, input);
				if (input >> line && line != "</V>") std::cout << "Some troubles in <V>" << std::endl;
			}
			else if (line == "<P>") {
				key_P = TRUE;
				Matrix_resize(P_n, par.N1_p, par.N2_p);
				Matrix_read  (P_n, par.N1_p, par.N2_p, input);
				if (input >> line && line != "</P>") std::cout << "Some troubles in <P>" << std::endl;
			}
			else if (line == "<Solids>")
				Read_Solids(input, Solids, Nodes, par);
		}
	}
	else {
		std::cout << "Error: No file! Starting Problem with Default Parameters" << std::endl;
		std::getchar();
	}
	
	if (key_U == FALSE) Matrix_resize(U_n, par.N1_u, par.N2_u);
	if (key_V == FALSE) Matrix_resize(V_n, par.N1_v, par.N2_v);
	if (key_P == FALSE) Matrix_resize(P_n, par.N1_p, par.N2_p);
	
	if (par.BC == Lamb_Oseen) {
		for (auto& it : Solids) {
			it.omega_n[3] = Lamb_Oseen_omega(it.r, par.Re, 0.0, par.Lamb_Oseen_r0);
			it.omega[3] = Lamb_Oseen_omega(it.r, par.Re, par.d_t, par.Lamb_Oseen_r0);
		}
	}

	if (key_U == FALSE && key_V == FALSE && key_P == FALSE)
	fill_exact(U_n, V_n, P_n, par, 0.0, par.d_t*0.5);             // Initial conditions for velocity and pressure

	std::cout << "Loading Data Finished" << std::endl;

}

void Post(fs::path WorkDir, Param &par, std::vector<Solid>& Solids, std::vector<Node> &Nodes, Matrix& U_n, Matrix& V_n, Matrix& P_n) {

	std::string file = WorkDir.string() + "/" + "input.txt";
	Load_Data(file, par, Solids, Nodes, U_n, V_n, P_n);

	file = "history_post_average10_dist5";
	history_init(WorkDir.string(), file, par.BC);

	int i = 1;
	while (Read_plt(par.WorkDir + "/" + "step" + std::to_string(i * 100) + ".plt", par, Solids)) {
		std::cout << "i = " << i * 100 << ", start new file" << std::endl;

		double h_average;
		h_average_of_Solids_Layer(Solids, par, h_average);
		history_log(par.WorkDir, file, par.N_step*par.d_t, h_average, 0., 0.);
		Solids.clear();
		i++;
	}
	std::cin.ignore();
}

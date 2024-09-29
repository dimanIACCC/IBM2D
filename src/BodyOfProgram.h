#pragma once
#include "stdafx.h"
#include "CalculateForce.h"
#include "Calculate_u_p.h"
#include "Output.h"
#include "PredictVel.h"


void BodyOfProgram(Param& par, std::vector<Solid>& solidList, std::vector<Node>& Nodes, Matrix& U_n, Matrix& V_n, Matrix& P);
void Save_Data(Param& par, std::vector<Solid>& solidList, std::vector<Node>& Nodes, Matrix& U_n, Matrix& V_n, Matrix& P);
void Load_Data(std::string &file, Param& par, std::vector<Solid>& solidList, std::vector<Node>& Nodes, Matrix& U_n, Matrix& V_n, Matrix& P);
void Post(fs::path WorkDir, Param &par, std::vector<Solid>& Solids, std::vector<Node> &Nodes, Matrix& U_n, Matrix& V_n, Matrix& P_n);

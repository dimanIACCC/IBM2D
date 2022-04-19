#pragma once
#include "stdafx.h"
#include "CalculateForce.h"
#include "Calculate_u_p.h"
#include "Output.h"
#include "PredictVel.h"


void BodyOfProgram(Param& par, std::vector<Circle>& solidList, std::vector<Node>& Nodes, Matrix& U_n, Matrix& V_n, Matrix& P);
void MakeHibernationFile(Param& par, std::vector<Circle>& solidList, std::vector<Node>& Nodes, Matrix& U_n, Matrix& V_n, Matrix& P);
void Awake(std::string &file, Param& par, std::vector<Circle>& solidList, std::vector<Node>& Nodes, Matrix& U_n, Matrix& V_n, Matrix& P);

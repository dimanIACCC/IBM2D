#include "PredictVel.h"


void ExplicPredVel(Matrix& U_predict, Matrix& V_predict, Matrix& U_n, Matrix& V_n, Matrix& P, Matrix& Force_x, Matrix& Force_y, Grid grid) {
	
	for (int i = 0; i < U_predict.size(); i++) {
		U_predict[i][0]                    = U_n[i][0]; 
		U_predict[i][U_predict[0].size() - 1] = U_n[i][U_predict[0].size() -1];
	}
	//for (int i = 0; i < U_predict[0].size(); i++) { U_predict[0][i] = U_n[0][i]; U_predict[U_predict[0].size() -1 ][i] = U_n[U_predict[0].size() -1][i];}
	for (int i = 0; i < V_predict.size(); i++) {
		V_predict[i][0] = V_n[i][0];
		V_predict[i][V_predict[0].size() - 1] = V_n[i][V_predict[0].size() - 1];
	}
	//for (int i = 0; i < U_predict[0].size(); i++);
	

	for (int i = 1; i < U_predict.size()-1; i++) {
		for (int j = 1; j < U_predict[0].size()-1; j++)
		{
			double LaplasU = (U_n[i + 1][j] - 2 * U_n[i][j] + U_n[i - 1][j]) / pow(grid.d_x, 2) + (U_n[i][j + 1] - 2 * U_n[i][j] + U_n[i][j - 1]) / pow(grid.d_y, 2);
			double CentrDiffx_U = (U_n[i + 1][j] - U_n[i - 1][j]) / (2 * grid.d_x);
			double CentrDiffy_U = (U_n[i][j + 1] - U_n[i][j - 1]) / (2 * grid.d_y);
			double GradPressX   = (P[i + 1][j] - P[i - 1][j]) / (2 * grid.d_x);
			U_predict[i][j] = grid.d_t*(LaplasU - U_n[i][j] * CentrDiffx_U - V_n[i][j] * CentrDiffy_U) + U_n[i][j] - GradPressX +Force_x[i][j];
		}
	}
	for (int i = 1; i < V_predict.size()-1; i++) {
		for (int j = 1; j < V_predict[0].size()-1; j++)
		{
			double LaplasV = (V_n[i + 1][j] - 2 * V_n[i][j] + V_n[i - 1][j]) / pow(grid.d_x, 2) + (V_n[i][j + 1] - 2 * V_n[i][j] + V_n[i][j - 1]) / pow(grid.d_y, 2);
			double CentrDiffx_V = (V_n[i + 1][j] - V_n[i - 1][j]) / (2 * grid.d_x);
			double CentrDiffy_V = (V_n[i][j + 1] - V_n[i][j - 1]) / (2 * grid.d_y);
			double GradPressY   = (P[i][j + 1] - P[i][j - 1]) / (2 * grid.d_y);
			V_predict[i][j] = grid.d_t*(LaplasV - U_n[i][j] * CentrDiffx_V - V_n[i][j] * CentrDiffy_V) + V_n[i][j] - GradPressY +Force_y[i][j];
		}
	}



}
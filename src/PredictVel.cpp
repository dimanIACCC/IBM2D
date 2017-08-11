#include "PredictVel.h"


void ExplicPredVel(Matrix& U_predict, Matrix& V_predict, Matrix& U_n, Matrix& V_n, int const n1, int const n2, Matrix operator_A[5], Matrix &b, Grid grid) {
	
	for (int i = 0; i < n1; i++) {
		for (int j = 0; j < n2; j++)
		{
			double LaplasU = (U_n[i + 1][j] - 2 * U_n[i][j] + U_n[i - 1][j]) / pow(grid.d_x, 2) + (U_n[i][j + 1] - 2 * U_n[i][j] + U_n[i][j - 1]) / pow(grid.d_y, 2);
			double LaplasV = (V_n[i + 1][j] - 2 * V_n[i][j] + V_n[i - 1][j]) / pow(grid.d_x, 2) + (V_n[i][j + 1] - 2 * V_n[i][j] + V_n[i][j - 1]) / pow(grid.d_y, 2);
			double ForwardDiffx_U = (U_n[i + 1][j] - U_n[i][j]) / grid.d_x;
			double ForwardDiffy_U = (U_n[i][j + 1] - U_n[i][j]) / grid.d_y;
			double ForwardDiffx_V = (V_n[i + 1][j] - V_n[i][j]) / grid.d_x;
			double ForwardDiffy_V = (V_n[i][j + 1] - V_n[i][j]) / grid.d_y;
			U_predict[i][j] = grid.d_t*(LaplasU - U_n[i][j] * ForwardDiffx_U - V_n[i][j] * ForwardDiffy_U) + U_n[i][j];
			V_predict[i][j] = grid.d_t*(LaplasV - U_n[i][j] * ForwardDiffx_V - V_n[i][j] * ForwardDiffy_V) + V_n[i][j];
		}
	}

	//for (int i = 0; i < n1; i++) {
	//	for (int j = 0; j < n2; j++)
	//	{
	//		double BackDiffx = (U_predict[i][j] - U_predict[i - 1][j]) / grid.d_x;
	//		double BackDiffy = (V_predict[i][j] - V_predict[i][j - 1]) / grid.d_y;

	//		U_predict[i][j] -= BackDiffx + BackDiffy;
	//		V_predict[i][j] = grid.d_t*(LaplasV - U_n[i][j] * ForwardDiffx_V - V_n[i][j] * ForwardDiffy_V) + V_n[i][j];
	//	}
	//}

}
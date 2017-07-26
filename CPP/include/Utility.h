#ifndef UTILITY_H
#define UTILITY_H

#include "mex.h"
#include "armadillo"
#include "gurobi_c++.h"

class Utility
{
protected:
	// MATLAB <--> Amadillo:
	static void ConvertSparseMatrixMatlabToArmadillo(const mxArray* mtlbMatrix, arma::sp_mat* armaMatrix);
	static void ConvertDenseMatrixMatlabToArmadillo(const mxArray* mtlbMatrix, arma::mat* armaMatrix);
	static void ConvertSparseMatrixArmadilloToMatlab(const arma::sp_mat* armaMatrix, mxArray* mlbMatrix);
	static void ConvertDenseMatrixArmadilloToMatlab(const arma::mat* armaMatrix, mxArray* mlbMatrix);

	// Sparse <--> dense
	static void ConvertSparseMatrixToDense(const arma::sp_mat* sparseMatrix, arma::mat* denseMatrix);
	static void ConvertSparseMatrixToDense(const mxArray* sparseMatrix, mxArray* denseMatrix);
	static void ConvertDenseMatrixToSparse(const mxArray* denseMatrix, mxArray* sparseMatrix);
	static void ConvertDenseMatrixToSparse(const arma::mat* denseMatrix, arma::sp_mat* sparseMatrix);

public:
	// MATLAB model to Gurobi model:
	static void ReadMxModel(
		const mxArray* mxModel,
		int *numMetsPtr,
		int *numRxnsPtr,
		std::vector<double> *lb,
		std::vector<double> *ub,
		arma::sp_mat *S);

	static void CreateGurobiEnv(
		GRBEnv *env,
		const mxArray* mxGurobiOptions);

	static void CreateGurobiModel(
		GRBModel *model,
		int numMets,
		int numRxns,
		double *lb,
		double *ub,
		arma::sp_mat *SPtr);

	// MATLAB cell arrays:
	//static void ReadMxCellArray(mxArray* cellArray, std::vector<std::string>* str);
	//static void WriteMxCellArray(std::vector<std::string>* str, mxArray* cellArray);

	// MATLAB <--> Amadillo:
	static void ConvertMatrix(const mxArray* mtbMatrixIn, arma::sp_mat* armaMatrixOut);
	static void ConvertMatrix(const mxArray* mtbMatrixIn, arma::mat* armaMatrixOut);
	static void ConvertMatrix(const arma::mat* armaMatrixIn, mxArray* mtbMatrixOut);
	static void ConvertMatrix(const arma::sp_mat* armaMatrixIn, mxArray* mtbMatrixOut);

	static void mxPrintVec(arma::vec* v);
	static void mxPrintMat(arma::mat* A);
};

#endif // UTILITY_H
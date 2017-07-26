#include "Utility.h"

void Utility::ReadMxModel(
	const mxArray* mxModel,
	int *numMetsPtr,
	int *numRxnsPtr,
	std::vector<double> *lb,
	std::vector<double> *ub,
	arma::sp_mat *S)
{
	// Extract model fields:
	mxArray* mxS = mxGetField(mxModel, 0, "S");
	mxArray* mxLb = mxGetField(mxModel, 0, "lb");
	mxArray* mxUb = mxGetField(mxModel, 0, "ub");

	int numMets = (int)mxGetM(mxS);
	int numRxns = (int)mxGetN(mxS);

	numMetsPtr[0] = numMets;
	numRxnsPtr[0] = numRxns;

	lb->resize(numRxns);
	ub->resize(numRxns);

	// Check array sizes:
	if (mxLb == NULL || mxGetM(mxLb) != numRxns || mxGetN(mxLb) != 1)
		mexErrMsgTxt("lb must be a column vector of length 'numRxns'.");
	if (mxUb == NULL || mxGetM(mxUb) != numRxns || mxGetN(mxUb) != 1)
		mexErrMsgTxt("ub must be a column vector of length 'numRxns'.");

	double *mxLbPtr = mxGetPr(mxLb); // Lower bounds
	double *mxUbPtr = mxGetPr(mxUb); // Upper bounds

	std::copy(mxLbPtr, mxLbPtr + numRxns, &*lb->begin());
	std::copy(mxUbPtr, mxUbPtr + numRxns, &*ub->begin());

	Utility::ConvertMatrix(mxS, S); // S
}

void Utility::CreateGurobiEnv(
	GRBEnv *env,
	const mxArray* mxGurobiOptions)
{
	// Extract fields:
	mxArray* mxFeasTol = mxGetField(mxGurobiOptions, 0, "feasibilityTol");
	mxArray* mxOptTol = mxGetField(mxGurobiOptions, 0, "optimalityTol");
	mxArray* mxTimeLimit = mxGetField(mxGurobiOptions, 0, "timeLimit");
	mxArray* mxScaleFlag = mxGetField(mxGurobiOptions, 0, "scaleFlag");
	mxArray* mxMethod = mxGetField(mxGurobiOptions, 0, "method");
	mxArray* mxOutputFlag = mxGetField(mxGurobiOptions, 0, "outputFlag");

	// Set options:
	if (mxFeasTol != NULL && mxGetM(mxFeasTol) != 0)
		env->set(GRB_DoubleParam_FeasibilityTol, ((double *)mxGetPr(mxFeasTol))[0]);
	if (mxOptTol != NULL && mxGetM(mxOptTol) != 0)
		env->set(GRB_DoubleParam_OptimalityTol, ((double *)mxGetPr(mxOptTol))[0]);
	if (mxTimeLimit != NULL && mxGetM(mxTimeLimit) != 0)
		env->set(GRB_DoubleParam_TimeLimit, ((double *)mxGetPr(mxTimeLimit))[0]);
	if (mxScaleFlag != NULL && mxGetM(mxScaleFlag) != 0)
		env->set(GRB_IntParam_ScaleFlag, ((int *)mxGetPr(mxScaleFlag))[0]);
	if (mxMethod != NULL && mxGetM(mxMethod) != 0)
		env->set(GRB_IntParam_Method, ((int *)mxGetPr(mxMethod))[0]);
	if (mxOutputFlag != NULL && mxGetM(mxOutputFlag) != 0)
		env->set(GRB_IntParam_OutputFlag, ((int *)mxGetPr(mxOutputFlag))[0]);

	// Print options:
	//mexPrintf("\nOptions:\n");
	//mexPrintf("FeasibilityTol = %f\n", env->get(GRB_DoubleParam_FeasibilityTol));
	//mexPrintf("OptimalityTol = %f\n", env->get(GRB_DoubleParam_OptimalityTol));
	//mexPrintf("TimeLimit = %f\n", env->get(GRB_DoubleParam_TimeLimit));
	//mexPrintf("ScaleFlag = %d\n", env->get(GRB_IntParam_ScaleFlag));
	//mexPrintf("Method = %d\n", env->get(GRB_IntParam_Method));
	//mexPrintf("OutputFlag = %d\n\n", env->get(GRB_IntParam_OutputFlag));
}

void Utility::CreateGurobiModel(
	GRBModel *model,
	int numMets,
	int numRxns,
	double *lb,
	double *ub,
	arma::sp_mat *S)
{
	model->addVars(lb, ub, NULL, NULL, NULL, numRxns);

	model->update(); // Integrate variables

	std::vector<GRBLinExpr> exprs(numMets); // Gurobi constraint expressions (/metabolites)

	GRBVar *vars = model->getVars();

	arma::sp_mat::const_iterator start = S->begin();
	arma::sp_mat::const_iterator end = S->end();

	for (arma::sp_mat::const_iterator it = start; it != end; ++it) // Iterate over sparse stoichiometry matrix
		exprs.at(it.row()) += *it * vars[it.col()];

	std::vector<char> senses(numMets, GRB_EQUAL);
	std::vector<double> rhs(numMets, 0);
	std::vector<std::string> names(numMets, "");

	model->addConstrs(&exprs[0], &senses[0], &rhs[0], NULL, numMets);

	model->update(); // Integrate constraints
}

/*
 * Converts a sparse or dense mxArray matrix to a sparse armadillo matrix
 */
void Utility::ConvertMatrix(const mxArray* mtbMatrixIn, arma::sp_mat* armaMatrixOut)
{
	if (!mxIsSparse(mtbMatrixIn))
	{
		mxArray* tmp = mxCreateSparse(0, 0, 0, mxREAL);
		Utility::ConvertDenseMatrixToSparse(mtbMatrixIn, tmp);
		Utility::ConvertSparseMatrixMatlabToArmadillo(tmp, armaMatrixOut);
		mxDestroyArray(tmp);
	}
	else
	{
		Utility::ConvertSparseMatrixMatlabToArmadillo(mtbMatrixIn, armaMatrixOut);
	}
}

/*
 * Converts a sparse or dense mxArray matrix to a dense armadillo matrix
 */
void Utility::ConvertMatrix(const mxArray* mtbMatrixIn, arma::mat* armaMatrixOut)
{
	if (mxIsSparse(mtbMatrixIn))
	{
		mxArray* tmp = mxCreateDoubleMatrix(0, 0, mxREAL);
		Utility::ConvertSparseMatrixToDense(mtbMatrixIn, tmp);
		Utility::ConvertDenseMatrixMatlabToArmadillo(tmp, armaMatrixOut);
		mxDestroyArray(tmp);
	}
	else
	{
		Utility::ConvertDenseMatrixMatlabToArmadillo(mtbMatrixIn, armaMatrixOut);
	}
}

/*
 * Converts an armadillo matrix to Matlab sparse or dense
 */
void Utility::ConvertMatrix(const arma::mat* armaMatrixIn, mxArray* mtbMatrixOut)
{
	if (mxIsSparse(mtbMatrixOut))
	{
		arma::sp_mat tmp; 
		Utility::ConvertDenseMatrixToSparse(armaMatrixIn, &tmp);
		Utility::ConvertSparseMatrixArmadilloToMatlab(&tmp, mtbMatrixOut);
	}
	else
	{
		Utility::ConvertDenseMatrixArmadilloToMatlab(armaMatrixIn, mtbMatrixOut);
	}
}

/*
 * Converts a sparse armadillo matrix to Matlab sparse or dense
 */
void Utility::ConvertMatrix(const arma::sp_mat* armaMatrixIn, mxArray* mtbMatrixOut)
{
	if (mxIsSparse(mtbMatrixOut))
	{
		Utility::ConvertSparseMatrixArmadilloToMatlab(armaMatrixIn, mtbMatrixOut);
	}
	else
	{
		arma::mat tmp; 
		Utility::ConvertSparseMatrixToDense(armaMatrixIn, &tmp);
		Utility::ConvertDenseMatrixArmadilloToMatlab(&tmp, mtbMatrixOut);
	}
}

/*
 * Converts a sparse Matlab matrix to a sparse Armadillo matrix 
 */
void Utility::ConvertSparseMatrixMatlabToArmadillo(const mxArray* mtlbMatrix, arma::sp_mat* armaMatrix)
{
	if (!mxIsSparse(mtlbMatrix))
	{
		mexErrMsgTxt("The supplied matrix is not in sparse format \n");
		return;
	}

	arma::mat tmp;
	Utility::ConvertMatrix(mtlbMatrix, &tmp);

	(*armaMatrix) = tmp;
}

/*
 * Converts a dense Matlab matrix to a dense Armadillo matrix 
 */
void Utility::ConvertDenseMatrixMatlabToArmadillo(const mxArray* mtlbMatrix, arma::mat* armaMatrix)
{
	// Things from Matlab (aka source):
	int nRows = (int)mxGetM(mtlbMatrix);
	int nCols = (int)mxGetN(mtlbMatrix);
	double* srcPtr = mxGetPr(mtlbMatrix);

	if (mxIsSparse(mtlbMatrix))
	{
		mexErrMsgTxt("The supplied matrix is sparse while it must be dense.\n");

		return;
	}

	// Resize arma matrix:
	armaMatrix->resize(nRows, nCols); 

	// Copy data:
	double* destPtr = armaMatrix->memptr();
	std::copy(srcPtr, srcPtr + nRows*nCols, destPtr);
}

/*
 * Converts a sparse Armadillo matrix to a sparse Matlab matrix 
 */
void Utility::ConvertSparseMatrixArmadilloToMatlab(const arma::sp_mat* armaMatrix, mxArray* mlbMatrix)
{
	const arma::uword nMets = armaMatrix->n_rows;
	const arma::uword nRxns = armaMatrix->n_cols;
	const arma::uword nnz = armaMatrix->n_nonzero;

	// Re-allocate memory: 
	mxSetM(mlbMatrix, nMets);
	mxSetN(mlbMatrix, nRxns);
	mxSetNzmax(mlbMatrix, nnz);
	
	// Allocate memory for row and column pointers:
	mwIndex* jc = (mwIndex*)mxCalloc(nRxns + 1, sizeof(mwIndex));
	mwIndex* ir = (mwIndex*)mxCalloc(nnz, sizeof(mwIndex));
	double* pr = (double*)mxCalloc(nnz, sizeof(double));
	
	// Copy values:
	std::copy(armaMatrix->col_ptrs, armaMatrix->col_ptrs + nRxns+1, jc);
	std::copy(armaMatrix->row_indices, armaMatrix->row_indices + nnz, ir);
	std::copy(armaMatrix->values, armaMatrix->values + nnz, pr);

	// Store sparse info in sparse mxArray:
	mxSetJc(mlbMatrix, jc);
	mxSetIr(mlbMatrix, ir);
	mxSetPr(mlbMatrix, pr);
}

/*
 * Converts a dense Armadillo matrix to a dense Matlab matrix 
 */
void Utility::ConvertDenseMatrixArmadilloToMatlab(const arma::mat* armaMatrix, mxArray* mlbMatrix)
{
	// Get pointers:
	const double* armaPtr = armaMatrix->memptr();
	double* mlbPtr = (double*)mxMalloc(armaMatrix->n_elem * sizeof(double));

	// Copy data:
	std::copy(armaPtr, armaPtr + armaMatrix->n_elem, mlbPtr);
	mxSetM(mlbMatrix, armaMatrix->n_rows);
	mxSetN(mlbMatrix, armaMatrix->n_cols);
	mxSetPr(mlbMatrix, mlbPtr);
}

/*
 * Converts a sparse armadillo matrix to a dense matrix
 */
void Utility::ConvertSparseMatrixToDense(const arma::sp_mat* sparseMatrix, arma::mat* denseMatrix)
{
	(*denseMatrix) = (*sparseMatrix);
}

/*
 * Converts a sparse matlab matrix to a dense matrix
 */
void Utility::ConvertSparseMatrixToDense(const mxArray* sparseMatrix, mxArray* denseMatrix)
{
	// Resize and fill with zeros:
	int nRows = (int)mxGetM(sparseMatrix);
	int nCols = (int)mxGetN(sparseMatrix);

	// Pointers to sparse ingoing matrix:
	mwIndex* jc = mxGetJc(sparseMatrix);
	mwIndex* ir = mxGetIr(sparseMatrix);
	double* prSparse = mxGetPr(sparseMatrix);

	// Allocate memory:
	double* prDense = (double*)mxCalloc(nRows * nCols, sizeof(double));

	mwIndex k = 0;

	for (mwIndex i = 0; i < (mwIndex)nCols; i++)
	{
		for (mwIndex j = jc[i]; j < jc[i + 1]; j++)
		{
			prDense[nRows*i + ir[j]] = prSparse[k]; // Replace non-zero entries

			k++;
		}
	}

	// Pointers to dense outgoing matrix:
	mxSetM(denseMatrix, nRows);
	mxSetN(denseMatrix, nCols);
	mxSetPr(denseMatrix, prDense);
}

/*
 * Converts a dense matlab matrix to a sparse matrix
 */
void Utility::ConvertDenseMatrixToSparse(const mxArray* denseMatrix, mxArray* sparseMatrix)
{
	// Get the matrix dimensions:
	int nRows = (int)mxGetM(denseMatrix);
	int nCols = (int)mxGetN(denseMatrix);

	// Get a pointer to the data:
	double* prIn = mxGetPr(denseMatrix);

	// Allocate space and temporary vectors:
	mwIndex* jcOut = (mwIndex*)mxCalloc(nCols+1, sizeof(mwIndex));
	std::vector<double> prTmp;
	std::vector<mwIndex> irTmp;
	
	mwIndex nnz = 0;

	// Get the non-zero values:
	for (mwIndex i=0; i < (mwIndex)nCols; i++)
	{
		for (mwIndex j=0; j < (mwIndex)nRows; j++)
		{
			mwIndex index = i*nRows + j;

			if (prIn[index] != 0)
			{
				// Add only non-zero elements:
				prTmp.push_back(prIn[index]); 
				irTmp.push_back(j);

				nnz++;
			}
		}

		// Update jc array:
		jcOut[i + 1] = nnz;
	}

	// Allocate memory:
	mwIndex* irOut = (mwIndex*)mxCalloc(nnz, sizeof(mwIndex));
	double* prOut = (double*)mxCalloc(nnz, sizeof(double) );
	
	// Copy values:
	std::copy(&irTmp[0], &irTmp[0] + nnz, irOut);
	std::copy(&prTmp[0], &prTmp[0] + nnz, prOut);

	// Set matrix dimensions and data:
	mxSetM(sparseMatrix, nRows);
	mxSetN(sparseMatrix, nCols);
	mxSetJc(sparseMatrix, jcOut);
	mxSetIr(sparseMatrix, irOut);
	mxSetPr(sparseMatrix, prOut);
}

/*
 * Converts a dense armadillo matrix to a sparse matrix
 */
void Utility::ConvertDenseMatrixToSparse(const arma::mat* denseMatrix, arma::sp_mat* sparseMatrix)
{
	(*sparseMatrix) = (*denseMatrix);
}

void Utility::mxPrintVec(arma::vec* v)
{
	double* vPtr = v->memptr();

	mexPrintf("\n");

	for (int i = 0; i < v->size(); i++)
	{
		double val = vPtr[i];

		if (val >= 0) {	mexPrintf("    %f\n", val); } else { mexPrintf("   %f\n", val); }
	}

	mexPrintf("\n");
}

void Utility::mxPrintMat(arma::mat* A)
{
	/*
	 * Note that the elements in matrix A are stored fortran-style in the memory:
	 * [ 0    m  ...  (n-1)*m ]
	 * [ 1   m+1 ... (n-1)*m+1]
	 * [ :    :          :    ]
	 * [m-1 2m-1 ...    n*m   ]
	 */

	double* APtr = A->memptr();

	int numRows = A->n_rows;
	int numCols = A->n_cols;

	int i;

	mexPrintf("\n");

	for (int iRow = 0; iRow < numRows; iRow++)
	{
		i = iRow;

		for (int iCol = 0; iCol < numCols; iCol++)
		{
			double val = APtr[i];

			if (val >= 0) {	mexPrintf("    %f", val); } else { mexPrintf("   %f", val); }

			i += numRows;
		}

		mexPrintf("\n");
	}

	mexPrintf("\n");
}
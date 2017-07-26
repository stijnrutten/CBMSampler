#include "Sampler.h"

#include "mex.h"
#include "blas.h"
#include "lapack.h"
#include "omp.h"
#include <time.h> // To seed the random number generator

#include "Utility.h"

#ifdef __cplusplus 
	extern "C" bool utIsInterruptPending();
#else
	extern bool utIsInterruptPending();
#endif

// Options:
const double ZERO_TOL = 1e-4; // Values below this are considered zero
const int NUM_START_TRIALS = 1000; // Maximum number of trials to start a chain
const int NUM_EXTEND_TRIALS = 1000; // Maximum number of trials to extend a chain
const int NUM_PROJECTIONS = 10; // Maximum number of trials to project x onto the null space

void Sampler(std::vector<arma::mat>* samplesPerThread, int numWarmups, arma::mat* warmups0, arma::mat* A, const arma::vec lb, const arma::vec ub, int numSteps, int numStepsBeforeProj, int numThreads)
{
	int numRxns = A->n_cols;

	/**********************
	 * Compute null space *
	 **********************/

	// SVD outputs:
	arma::mat U, Vt;
	arma::vec S;

	// Singular value decomposition:
	svd(&U, &S, &Vt, A);

	int rank = ((arma::uvec)find(S > ZERO_TOL)).size(); // Number of non-zero singular values

	// Nullspace; the columns of V corresponding to the zero singular values:
	arma::mat Nt = Vt.rows(rank, numRxns - 1);
	arma::mat N = Nt.t();

	arma::vec tmp(numRxns - rank, 1); // Used for projection onto null space (numRxns - rank = N.n_cols)

	/****************************************
	 * Move warmup samples towards centroid *
	 ****************************************/

	arma::mat warmups = *warmups0; // Copy inputted warmup samples

	double* warmupsPtr = warmups.memptr();

	arma::vec centroid = mean(warmups, 1); // Get centroid

	// Move warmup samples towards centroid:
	warmups = (1.0/3.0) * warmups + (2.0/3.0) * centroid * arma::ones(1, numWarmups); // Weighted average

	/********************
	 * Perform sampling *
	 ********************/

	srand(time(NULL)); // Initialize random seed

	std::vector<int> numDiscardedPerThread(numThreads, 0);

	double startTime = omp_get_wtime();

#pragma omp parallel num_threads(numThreads) default(none) firstprivate(numSteps, numStepsBeforeProj, numRxns, numWarmups, Nt, N, tmp, centroid) shared(samplesPerThread, warmupsPtr, lb, ub, numDiscardedPerThread)
	{
		int iThread = omp_get_thread_num(); // Get thread ID

		arma::mat *samplesCur = &samplesPerThread->at(iThread); // Current sample matrix

		int totalSamplesCur = samplesCur->n_cols; // Current total number of samples (including warmup samples)
		int numSamplesCur = numWarmups;

		std::copy(warmupsPtr, warmupsPtr + numRxns*numWarmups, samplesCur->memptr()); // Copy warmup samples

		int iExtendTrial = 1;

		while (numSamplesCur < totalSamplesCur)
		{
			/*********************************
			 * Check for MATLAB interruption *
			 *********************************/

			if (utIsInterruptPending())
				break;

			/*********************
			 * Get random sample *
			 *********************/

			arma::vec x = samplesCur->col(rand() % numSamplesCur); // Pick random sample

			/****************
			 * Extend chain *
			 ****************/

			bool valid = true;

			int numStepsCur = 1;
			int iStartTrial = 1;

			while (numStepsCur < numSteps)
			{
				/************************
				 * Get random direction *
				 ************************/

				arma::vec x0 = samplesCur->col(rand() % numSamplesCur); // Pick random sample

				arma::vec u = x0 - centroid;

				// Make sure the random sample is different from the centroid:
				if (arma::all(u == 0))
				{
					if (iStartTrial == NUM_START_TRIALS)
					{
						mexPrintf("Reached maximum number of trials permitted to start a new chain.\n");

						// Exit while loops:
						numStepsCur = numSteps;
						numSamplesCur = totalSamplesCur;
						valid = false;
					}

					iStartTrial++;
					numDiscardedPerThread[iThread]++;

					continue;
				}

				/********************
				 * Take random step *
				 ********************/

				arma::vec steps = arma::join_cols((ub - x) / u, (lb - x) / u);

				arma::uvec isNeg = arma::find(steps < 0);
				arma::uvec isPos = arma::find(steps > 0);

				double minStep, maxStep;

				// Minimal absolute negative step:
				if (isNeg.size() != 0)
					minStep = arma::max(steps(isNeg));
				else
					minStep = 0;

				// Minimal positive step:
				if (isPos.size() != 0)
					maxStep = arma::min(steps(isPos));
				else
					maxStep = 0;

				double step = ((double)rand() / (double)RAND_MAX) * (maxStep - minStep) + minStep; // Random step magnitude

				x += step * u; // Perform step

				numStepsCur++;

				/***************************
				 * Project onto null space *
				 ***************************/

				if (numStepsCur % numStepsBeforeProj == 0 || numStepsCur == numSteps)
				{
					// Project x onto null space: x = N * N^T * x
					MatrixMultiply(&tmp, &Nt, &x);
					MatrixMultiply(&x, &N, &tmp);

					if (arma::any(x < (lb - ZERO_TOL)) || arma::any(x >(ub + ZERO_TOL)))
					{
						valid = false;
						break; // Stop the chain
					}
				}
			}

			/*******************
			 * Validate sample *
			 *******************/

			if (!valid)
			{
				mexPrintf("Breaking chain at step %d.\n", numStepsCur);

				if (iExtendTrial == NUM_EXTEND_TRIALS)
				{
					mexPrintf("Reached maximum number of trials permitted to extend a chain.\n");

					numSamplesCur = totalSamplesCur; // Exit while loop
				}

				iExtendTrial++;
				numDiscardedPerThread[iThread]++;

				continue;
			}

			/******************
			 * Add new sample *
			 ******************/

			samplesCur->col(numSamplesCur) = x; // Add sample

			centroid = (numSamplesCur * centroid + x) / (numSamplesCur + 1); // Update centroid

			numSamplesCur++;

			iExtendTrial = 1;
		}
	}

	double endTime = omp_get_wtime();

	int numDiscarded = 0;

	for (int iThread = 0; iThread < numThreads; iThread++)
	{
		numDiscarded += numDiscardedPerThread[iThread];
	}

	/********************
	 * Sampling summary *
	 ********************/

	//mexPrintf(" ========================= CBM SAMPLER SUMMARY =========================\n");
	//mexPrintf("Sampling phase took %4.2f seconds.\n", endTime - startTime);
	//mexPrintf("%d samples outside the nullspace have been discarded and were resampled.\n", numDiscarded);
	//mexPrintf(" =======================================================================\n");
}

/*--------------------------------------------------------------*
 * Singular value decomposition by LAPACK                       *
 *--------------------------------------------------------------*/
void svd(arma::mat* U, arma::vec* S, arma::mat* Vt, arma::mat* A)
{
	ptrdiff_t numRows = A->n_rows;
	ptrdiff_t numCols = A->n_cols;
	ptrdiff_t LWORK = std::min(numRows, numCols) * (6 + 4 * std::min(numRows, numCols)) + std::max(numRows, numCols);

	// SVD outputs:
	S->resize(std::min((int)numRows, (int)numCols));
	U->resize((int)numRows, (int)numRows);
	Vt->resize((int)numCols, (int)numCols);

	std::vector<ptrdiff_t> IWORK(8 * std::min(numRows, numCols));
	ptrdiff_t INFO;
	std::vector<double> WORK(LWORK);

	// Lapack Singular Value Decomposition (SVD):
	dgesdd("A", &numRows, &numCols, A->memptr(), &numRows, S->memptr(), U->memptr(), &numRows, Vt->memptr(), &numCols, &WORK[0], &LWORK, &IWORK[0], &INFO);
}

/*--------------------------------------------------------------*
 * Matrix multiplication by BLAS                                *
 *--------------------------------------------------------------*/
void MatrixMultiply(arma::vec* result, arma::mat* A, arma::vec* b)
{
	// Options:
	double alpha = 1.0;
	double beta = 0.0;

	// Dimensions of input matrices:
	ptrdiff_t m = A->n_rows;
	ptrdiff_t k = A->n_cols;
	ptrdiff_t n = b->n_cols;

	char transp = 'N'; // Do not transpose

	// Perform matrix multiplication using DGEMM from ACML:
	dgemm(&transp, &transp, &m, &n, &k, &alpha, A->memptr(), &m, b->memptr(), &k, &beta, result->memptr(), &m);
}
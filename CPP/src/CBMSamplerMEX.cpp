#include "Utility.h"
#include "SamplerFVA.h"
#include "Sampler.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
/* MEX function for CBMSamplerCPP. Samples the null space of a given CB model.
 *
 * Input:
 * 1. mxModel: CB model.
 * 2. numSamples: Number of samples to generate.
 * 3. numSteps: Number of steps between samples.
 * 4. numStepsBeforeProj: Number of steps before projection.
 * 5. numThreads: Number of parallel threads.
 * 6. mxWarmups (OPTIONAL): Warm-up samples.
 * 7. mxGurobiOptions (OPTIONAL): Gurobi options
 *
 * Output:
 * 1. samples: generated samples.
 * 2. warmups: Used warm-up samples.
 */
{
	/*--------------------------------------------------------------*
	 * Input                                                        *
	 *--------------------------------------------------------------*/

	if (nrhs < 5)
		mexErrMsgTxt("Invalid number of input arguments.");

	const mxArray *mxModel = prhs[0]; // CB Model
	
	int numSamples = ((int *)mxGetData(prhs[1]))[0]; // Number of samples
	int numSteps = ((int *)mxGetData(prhs[2]))[0]; // Number of steps between each sample
	int numStepsBeforeProj = ((int *)mxGetData(prhs[3]))[0]; // Number of steps before projection into null space
	int numThreads = ((int *)mxGetData(prhs[4]))[0]; // Number of threads
	
	const mxArray *mxWarmups; // Warm-up samples
	if (nrhs >= 6)
		mxWarmups = prhs[5];

	const mxArray *mxGurobiOptions; // Gurobi options
	if (nrhs >= 7)
		mxGurobiOptions = prhs[6];

	bool doWarmup = (mxWarmups == NULL || mxGetM(mxWarmups) == 0);

	/*--------------------------------------------------------------*
	 * Read CB model                                                *
	 *--------------------------------------------------------------*/

	int numMets /* Number of metabolites (/constraints) */,
		numRxns /* Number of reactions (/variables) */;
	std::vector<double> lb /* Lower bounds */, ub /* Upper bound */;
	arma::sp_mat S /* Stoichiometry matrix */;

	Utility::ReadMxModel(mxModel, &numMets, &numRxns, &lb, &ub, &S);

	/*--------------------------------------------------------------*
	 * Get warm-up samples                                          *
	 *--------------------------------------------------------------*/

	int numWarmups;

	arma::mat warmups;
	double* warmupsPtr;

	if (doWarmup) // Generate warm-up samples
	{
		numWarmups = 2 * numRxns;

		warmups = arma::mat(numRxns, numWarmups, arma::fill::zeros);
		warmupsPtr = warmups.memptr();

		// Create Gurobi environment:
		GRBEnv env = GRBEnv();
		Utility::CreateGurobiEnv(&env, mxGurobiOptions);

		// Create Gurobi model:
		GRBModel model = GRBModel(env);
		Utility::CreateGurobiModel(&model, numMets, numRxns, &lb[0], &ub[0], &S);

		SamplerFVA( // Perform FVA
			&model,
			&warmups,
			numThreads);
	}
	else // Read warm-up samples
	{
		numWarmups = (int)(mxGetDimensions(mxWarmups)[1]);

		double* mxWarmupsInPtr = (double *)mxGetData(mxWarmups);

		warmups = arma::mat(numRxns, numWarmups, arma::fill::zeros);
		warmupsPtr = warmups.memptr();

		std::copy(mxWarmupsInPtr, mxWarmupsInPtr + numRxns*numWarmups, warmupsPtr);
	}

	/*--------------------------------------------------------------*
	 * Perform sampling                                             *
	 *--------------------------------------------------------------*/

	std::vector<arma::mat> samplesPerThread(numThreads); // Vector containing the sample matrices

	if (numSamples > 0)
	{
		double numSamplesPerThread = (double)numSamples / (double)numThreads;

		int total = 0; // Total number of assigned samples

		for (int iThread = 0; iThread < numThreads; iThread++)
		{
			int numSamplesCur = (int)round((iThread + 1) * numSamplesPerThread) - total;

			total += numSamplesCur;

			arma::mat samplesCur(numRxns, numWarmups + numSamplesCur, arma::fill::zeros); // Contains warmup samples

			samplesPerThread[iThread] = samplesCur;
		}

		arma::mat A(S); // Turn sparse matrix into full matrix

		arma::vec lbArma(lb);
		arma::vec ubArma(ub);

		Sampler( // Perform sampling
			&samplesPerThread,
			numWarmups,
			&warmups,
			&A,
			lbArma,
			ubArma,
			numSteps,
			numStepsBeforeProj,
			numThreads);
	}

	/*--------------------------------------------------------------*
	 * Write output                                                 *
	 *--------------------------------------------------------------*/

	if (nlhs >= 1) // Samples
	{
		mxArray* mxSamples;

		if (numSamples > 0)
		{
			mxSamples = mxCreateDoubleMatrix(numRxns, numSamples, mxREAL);
			double* mxSamplesPtr = mxGetPr(mxSamples);

			int numStored = 0; // Number of sampled stored

			for (int iThread = 0; iThread < numThreads; iThread++)
			{
				arma::mat samplesCur = samplesPerThread[iThread];

				double* samplesCurPtr = samplesCur.memptr();

				int numSamplesCur = samplesCur.n_cols - numWarmups;

				std::copy(samplesCurPtr + numRxns*numWarmups, samplesCurPtr + numRxns*(numWarmups + numSamplesCur), mxSamplesPtr + numRxns*numStored); // Only copy samples (not warmups)

				numStored += numSamplesCur;
			}
		}
		else
		{
			mxSamples = mxCreateDoubleMatrix(0, 0, mxREAL); // Empty MEX array
		}

		plhs[0] = mxSamples;
	}

	if (nlhs >= 2) // Warmups
	{
		mxArray* mxWarmups = mxCreateDoubleMatrix(numRxns, numWarmups, mxREAL);
		double* mxWarmupsPtr = mxGetPr(mxWarmups);

		std::copy(warmupsPtr, warmupsPtr + numRxns*numWarmups, mxWarmupsPtr);

		plhs[1] = mxWarmups;
	}
}
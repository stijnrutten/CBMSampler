#include "SamplerFVA.h"

#include "omp.h"

void SamplerFVA(GRBModel* model, arma::mat* warmups, int numThreads)
{
	int numRxns = model->get(GRB_IntAttr_NumVars);

	arma::vec sample(numRxns);

	GRBEnv env = model->getEnv();
	GRBVar* vars = model->getVars();

	env.set(GRB_IntParam_Presolve, 2); // Aggressive presolve

	model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE); // Minimize fluxes

	GRBVar var0 = vars[0];

	var0.set(GRB_DoubleAttr_Obj, 1.0);
	model->optimize(); // Minimize model
	var0.set(GRB_DoubleAttr_Obj, 0.0);

	double* xVals = model->get(GRB_DoubleAttr_X, vars, numRxns);

	std::copy(
		xVals,
		xVals + numRxns,
		sample.memptr());

	warmups->col(0) = sample;

	env.set(GRB_IntParam_Presolve, 0); // Disable presolve

	if (numThreads == 1) // Perform FVA in a single thread
	{
		int iSample = 1;

		for (int iRxn = 1; iRxn < numRxns; iRxn++) // For each reaction
		{
			GRBVar var = vars[iRxn]; // Current reaction

			var.set(GRB_DoubleAttr_Obj, 1.0);
			
			model->optimize(); // Minimize model

			if (model->get(GRB_IntAttr_Status) != 2) // Not optimal: try again without warm basis
			{
				model->reset(); // Reset model to unsolved state
				model->optimize(); // Optimize model
			}

			if (model->get(GRB_IntAttr_Status) == 2) // Optimal
			{
				xVals = model->get(GRB_DoubleAttr_X, vars, numRxns);
				std::copy(xVals, xVals + numRxns, sample.memptr());
				warmups->col(iSample) = sample;
			}
			else // Failed
			{
			}

			var.set(GRB_DoubleAttr_Obj, 0.0);

			iSample++;
		}

		model->set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE); // Maximize fluxes

		for (int iRxn = 0; iRxn < numRxns; iRxn++)
		{
			GRBVar var = vars[iRxn];

			var.set(GRB_DoubleAttr_Obj, 1.0);

			model->optimize(); // Maximize model

			if (model->get(GRB_IntAttr_Status) != 2) // Not optimal: try again without warm basis
			{
				model->reset(); // Reset model to unsolved state
				model->optimize(); // Optimize model
			}

			if (model->get(GRB_IntAttr_Status) == 2) // Optimal
			{
				xVals = model->get(GRB_DoubleAttr_X, vars, numRxns);
				std::copy(xVals, xVals + numRxns, sample.memptr());
				warmups->col(iSample) = sample;
			}
			else // Failed
			{
			}

			var.set(GRB_DoubleAttr_Obj, 0.0);

			iSample++;
		}
	}
	else // Perform FVA in multiple threads
	{
		/*--------------------------------------------------------------*
		 * Create a copy of the model for each thread                   *
		 *--------------------------------------------------------------*/

		model->update(); // Update objective before copying the model

		// Models per thread:
		std::vector<GRBModel> models;
		models.reserve(numThreads);

		models.push_back(*model); // First model

		for (int iCopy = 1; iCopy < numThreads; iCopy++) // For each thread
			models.push_back(GRBModel(*model)); // Copy model

		/*--------------------------------------------------------------*
		 * Start parallel pool                                          *
		 *--------------------------------------------------------------*/

		if (numThreads % 2 == 0) // Even number of threads
		{
#pragma omp parallel num_threads(numThreads) default(none) private(xVals) firstprivate(numRxns, numThreads) shared(models, warmups)
			{
				arma::vec sampleCur(numRxns);

				int iThread = omp_get_thread_num();

				GRBModel modelCur = models[iThread];
				GRBVar *vars = modelCur.getVars();

				bool minimize = iThread < (numThreads / 2);

				int sampleOffset, iPart;

				if (minimize)
				{
					sampleOffset = 0;
					iPart = iThread;
				}
				else
				{
					sampleOffset = numRxns;
					iPart = iThread - numThreads / 2;

					modelCur.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE); // Maximize fluxes
				}

				double partSize = (double)numRxns / ((double)numThreads / 2.0);

				int iRxn0 = (int)floor(partSize * (double)iPart);

				if (minimize && iRxn0 == 0)
					iRxn0 = 1;

				int iRxnEnd = (int)floor(partSize * (double)(iPart + 1));

				for (int iRxn = iRxn0; iRxn < iRxnEnd; iRxn++)
				{
					GRBVar var = vars[iRxn]; // Current variable

					var.set(GRB_DoubleAttr_Obj, 1.0); // Set objective

					modelCur.optimize(); // Optimize model

					if (modelCur.get(GRB_IntAttr_Status) != 2) // Not optimal: try again without warm basis
					{
						modelCur.reset(); // Reset model to unsolved state
						modelCur.optimize(); // Optimize model
					}

					if (modelCur.get(GRB_IntAttr_Status) == 2) // Optimal
					{
						xVals = modelCur.get(GRB_DoubleAttr_X, vars, numRxns);

						std::copy(
							xVals,
							xVals + numRxns,
							sampleCur.memptr());

						warmups->col(iRxn + sampleOffset) = sampleCur;
					}
					else // Failed
					{
					}

					var.set(GRB_DoubleAttr_Obj, 0.0); // Reset objective
				}
			}
		}
		else // Uneven number of threads > 1
		{
#pragma omp parallel num_threads(numThreads) default(none) private(xVals) firstprivate(numRxns) shared(models, warmups)
			{
				arma::vec sampleCur(numRxns);

				int iThread = omp_get_thread_num();

				GRBModel modelCur = models[iThread];
				GRBVar *vars = modelCur.getVars();

#pragma omp for
				for (int iSample = 1; iSample < 2 * numRxns; iSample++)
				{
					bool minimize = iSample < numRxns;

					int iRxn;

					if (minimize)
					{
						iRxn = iSample;

						modelCur.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE); // Minimize objective
					}
					else
					{
						iRxn = iSample - numRxns;

						modelCur.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE); // Maximize objective
					}

					GRBVar var = vars[iRxn];

					var.set(GRB_DoubleAttr_Obj, 1.0);

					modelCur.optimize(); // Optimize model

					if (modelCur.get(GRB_IntAttr_Status) != 2) // Not optimal: try again without warm basis
					{
						modelCur.reset(); // Reset model to unsolved state
						modelCur.optimize(); // Optimize model
					}

					if (modelCur.get(GRB_IntAttr_Status) == 2) // Optimal
					{
						xVals = modelCur.get(GRB_DoubleAttr_X, vars, numRxns);
						
						std::copy(
							xVals,
							xVals + numRxns,
							sampleCur.memptr());
						
						warmups->col(iSample) = sampleCur;
					}
					else // Failed
					{
					}

					var.set(GRB_DoubleAttr_Obj, 0.0);
				}
			}
		}
	}
}
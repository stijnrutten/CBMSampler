#ifndef SAMPLERFVA_H
#define SAMPLERFVA_H

#include "armadillo"
#include "gurobi_c++.h"

void SamplerFVA(
	GRBModel* model,
	arma::mat* warmups,
	int numThreads);

#endif // SAMPLERFVA_H
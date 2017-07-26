#ifndef SAMPLER_H
#define SAMPLER_H

#include "armadillo"

void Sampler(
	std::vector<arma::mat>* samplesPerThread,
	int numWarmups,
	arma::mat* warmups0,
	arma::mat* A,
	const arma::vec lb,
	const arma::vec ub,
	int numSteps,
	int numStepsBeforeProj,
	int numThreads);

void svd(
	arma::mat* U,
	arma::vec* S,
	arma::mat* Vt,
	arma::mat* A);

void MatrixMultiply(
	arma::vec* result,
	arma::mat* A,
	arma::vec* b);

#endif // SAMPLER_H
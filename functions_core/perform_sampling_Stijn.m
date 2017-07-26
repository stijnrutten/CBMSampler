function [sModel,string] = perform_sampling_Stijn(model,nSamples,nSteps,numStepsBeforeProj,nThreads,warmups,gurobiOptions)
% Perform sampling using CBMSampler

fprintf('\nPerforming sampling...\n');

tic;
[sModel.points, sModel.warmups] = CBMSamplerCPP(model, nSamples, nSteps, ...
    numStepsBeforeProj, nThreads, warmups, gurobiOptions);

TSampling = toc;

fprintf('Done! It took %f seconds.\n', TSampling);

% Generate name of current calculation
string = sprintf('nSteps%inSamples%i',nSteps,nSamples);
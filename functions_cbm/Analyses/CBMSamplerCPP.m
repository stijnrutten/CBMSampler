function [samples, warmups] = CBMSamplerCPP(model, numSamples, numSteps, numStepsBeforeProj, numThreads, warmups, gurobiOptions)

%% Input

if nargin < 3 || isempty(numSteps), numSteps = 10; end % Number of steps between each sample
if nargin < 4 || isempty(numStepsBeforeProj), numStepsBeforeProj = 10; end % Number of steps before projection into null space
if nargin < 5 || isempty(numThreads), numThreads = 1; end % Number of parallel threads
if nargin < 6, warmups = []; end % Warmup samples
if nargin < 7 || isempty(gurobiOptions), gurobiOptions = struct(); end % Gurobi options

%% Initialize

[numMets, numRxns] = size(model.S);

if ~isfield(model, 'b'), model.b = zeros(numMets, 1); end
if ~isfield(model, 'c') || any(model.c), model.c = zeros(numRxns, 1); end

%% Call MEX

[samples, warmups] = CBMSamplerMEX(model, int32(numSamples), int32(numSteps), ...
    int32(numStepsBeforeProj), int32(numThreads), warmups, gurobiOptions);

end
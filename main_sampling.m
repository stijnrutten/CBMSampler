%% main code of sampler convergence

clear all
clc

%% Set additional local library paths and define some parameters

[~, hostName] = system('hostname');

switch strtrim(hostName)
    case 'Z50'
        modelPath = strcat(pwd,'\metabolic_models');
        functionsCBMPath = strcat(pwd,'\functions_cbm');
        functionsCorePath = strcat(pwd,'\functions_core');
        functionsMiscPath = strcat(pwd,'\functions_misc');
        MEXPath = strcat(pwd,'\MEX');
        nThreads = 4;
    case 'BMT-PC02'
        modelPath = strcat(pwd,'\metabolic_models');
        functionsCBMPath = strcat(pwd,'\functions_cbm');
        functionsCorePath = strcat(pwd,'\functions_core');
        functionsMiscPath = strcat(pwd,'\functions_misc');
        MEXPath = strcat(pwd,'\MEX');
        nThreads = 2;
    case 'ToshibaPortege'
        modelPath = strcat(pwd,'\metabolic_models');
        functionsCBMPath = strcat(pwd,'\functions_cbm');
        functionsCorePath = strcat(pwd,'\functions_core');
        functionsMiscPath = strcat(pwd,'\functions_misc');
        MEXPath = strcat(pwd,'\MEX');
        nThreads = 2;
    otherwise
        warning('Unrecognized host name. Local files are not available.');
end

addpath(modelPath)
addpath(functionsCorePath)
addpath(functionsMiscPath)
addpath(genpath(functionsCBMPath))
addpath(MEXPath)

%% Load model

[model] = perform_load_model_Stijn(modelPath);

%% Add rxnBoundary vector to model
model.rxnBoundary = findExcRxns(model);

%% Sampler options

numStepsBeforeProj = 25;
gurobiOptions = GurobiOptions('ScaleFlag', 0, 'Method', 0, 'OutputFlag', 0);

warmups = [];

nSamples = 10;
nSteps = 25;

%% Perform sampling

% Perform sampling
[sModel,string] = perform_sampling_Stijn(model,nSamples,nSteps,numStepsBeforeProj,nThreads,warmups,gurobiOptions);

% Fix thermodynamically infeasible loops
[sModel] = perform_FixLoops_Stijn(model,sModel);
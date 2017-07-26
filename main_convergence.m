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

no_sampler = 0; % 1 if samples are loaded, 0 if samples have to be generated

%% Add rxnBoundary vector to model
model.rxnBoundary = findExcRxns(model);

%% Sampler options

numStepsBeforeProj = 25;
gurobiOptions = GurobiOptions('ScaleFlag', 0, 'Method', 0, 'OutputFlag', 0);

warmups = [];
nSamples = 10;
nSteps = 25;

%% Convergence options

c_threshold = 0.9; % convergence threshold
accuracy_range = 0.03; % desired accuracy range around the found most probable flux value for CSM convergence

% select convergence method
[convergence_method] = perform_select_convergence_method_Stijn();

%% Perform solution space convergence

% initialize sampling and convergence parameters
converged = zeros(size(model.rxns,1),1); % variable that indicates convergence for each reaction (zero = not converged, 1 = converged)
i = 0; % loop counter
c = 0; % convergence tracker indicator, 1 if convegrence is achieved for a sufficient number of reactions

while c == 0
    
    i = i+1;
    
    if no_sampler == 1
        string = x_string_names{i};
    else
        
        % Allocate warmup points from the first calculation to be used for upcoming calculations
        if i == 2
            warmups = sModel.warmups;
        end
        
        % Perform sampling
        [sModel,string] = perform_sampling_Stijn(model,nSamples,nSteps,numStepsBeforeProj,nThreads,warmups,gurobiOptions);
        
        % Fix thermodynamically infeasible loops
        [sModel] = perform_FixLoops_Stijn(model,sModel);
        
        % Rename data from sampling function
        x.(string) = sModel.points_fixed;
        x_string_names = fieldnames(x);
        
    end
    
    % select chosen convergence method
    if convergence_method == 1
        % Perform convergence
        [converged,c] = perform_convergence_STM_Stijn(string,x,i,converged,nSteps,nSamples,model,c_threshold,c);
    elseif convergence_method == 2
        if i == 1
            x_max = 0;
        end
        [converged,c,x_max] = perform_convergence_CSM_Stijn(accuracy_range,string,model,x,x_string_names,i,converged,c_threshold,c,x_max,nSteps,nSamples);
    elseif convergence_method == 3
        fprintf('No valid convergence method selected. Convergence analysis ended.\n');
        break
    end
    
    % Increase the amount of samples for the next loop
    [nSamples,nSteps] = perform_increase_parameters_Stijn(c,nSamples,nSteps);
    
end
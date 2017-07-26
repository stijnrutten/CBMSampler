function [X, fixed] = FixLoops(model, X)

%% Initialize

numSamples = size(X, 2);

fixed = false(numSamples, 1);

%% Fix thermodynamically infeasible loops

for iSample = 1:numSamples % For each sample
    x = X(:, iSample); % Current sample
    
    % Setup cycle free flux LP problem:
    [LPProblem, rxnInds] = SetupMinAbsFluxLPGurobi(model, x);

    % Solve cycle free flux LP problem:
    sol = SolveLPGurobi(LPProblem);
    
    if strcmp(sol.status, 'OPTIMAL')
        X(rxnInds, iSample) = sol.x; % Set fixed fluxes
        
        fixed(iSample) = true;
    end
end

%% Skipped samples

if any(~fixed)
    indsStr = sprintf(' %d', find(~fixed));
    
    warning('Non-optimal solution for samples%s. Skipped.', indsStr);
end

end
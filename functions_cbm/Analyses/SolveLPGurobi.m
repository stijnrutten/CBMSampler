function sol = SolveLPGurobi(LPProblem)
% LPProblem: struct, containing the fields:
%    'A': stoichiometry (numMets*numRxns): model.S
%    'b': model.b
%    'c': objective (numRxns*1)
%    'lb': lower bounds (numRxns*1)
%    'ub': upper bounds (numRxns*1)
%    'csense': constraint sense (numMets*1): '=' (equal), '>' (greater) or '<' (lesser).
%    'osense': objective sense: 'min' or 'max'.
%    'cbasis' (OPTIONAL):
%    'vbasis' (OPTIONAL):

params = struct(...
    'OutputFlag',      0, ... % No ouput
    'DisplayInterval', 1, ...
    'FeasibilityTol',  1e-9, ...
    'OptimalityTol',   1e-9, ...
    'Method',          -1); % Automatically choose algorithm

% Solve LP using gurobi MEX:
sol = gurobi(LPProblem, params);

end
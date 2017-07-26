function [LPProblem, rxnInds] = SetupMinAbsFluxLPGurobi(model, v0, rxnConst)
% Sets up LP problem of the form:
%        min: sum(|v_i|)
%              i
% subject to: S*v = 0                                  (c. 1)
%             0 <= v_i <= v_i^0 for i with v_i^0 >= 0  (c. 2)
%             v_i^0 <= v_i <= 0 for i with v_i^0 < 0   (c. 3)
%             v_j = v_j^0 for all exchange fluxes v_j  (c. 4)
%             lb <= v <= ub                            (c. 5)
% 
% Objective:
%     The sum over all absolute fluxes equals the sum over positive fluxes
% minus the sum over negative fluxes:
% sum(|v_i|) =        sum       (v_i) -       sum       (v_i)
%  i           i with v_i^0 >= 0        i with v_i^0 < 0
% 
% Constraints:
% c. 1: Equillibrium.
% c. 2, 3: Absolute fluxes cannot increase (only decrease). New bounds.
% c. 4: Exchange fluxes cannot change. Turn vars into constants.
% c. 5: Bound constraints.

%% Input

if nargin < 3, rxnConst = model.rxnBoundary; end

%% Options

boundTol = 1e-7;

%% Set exchange reeaction fluxes constant

% Indices of nonzero-flux reactions:
rxnInds = find(abs(v0) > 0);

% Select nonzero-flux reactions:
model = ExtractSubModel(model, rxnInds);

v0 = v0(rxnInds);
rxnConst = rxnConst(rxnInds);

numMets = size(model.S, 1);

%% Linear equations

LPProblem.A = model.S; % A

%% Objective

v0Pos = bsxfun(@and, ... % Positive-flux for non-constant reactions
    v0 >= 0, ... % Positive flux
    ~rxnConst); % Not constant

v0Neg = bsxfun(@and, ... % Negative-flux for non-constant reactions
    v0 < 0, ... % Negative flux
    ~rxnConst); % Not constant

LPProblem.obj = double(v0Pos - v0Neg); % Objective
LPProblem.modelsense = 'min'; % osense

%% Constraints

LPProblem.rhs = zeros(numMets, 1); % Right hand side
LPProblem.sense = repmat('=', numMets, 1); % csense

%% Bounds
% Assume: lb <= 0 when ub > 0 and
%         ub >= 0 when lb < 0
%         for internal reactions.

LPProblem.lb = v0 - boundTol; % Lower bounds
LPProblem.ub = v0 + boundTol; % Upper bounds

% Allow decrease of absolute flux for internal reactions:
LPProblem.lb(v0Pos) = 0;
LPProblem.ub(v0Neg) = 0;

end
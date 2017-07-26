function rxnFields = GetRxnFields(model)

%% Initialize

[numMets, numRxns] = size(model.S);

fieldNames = fieldnames(model);

%% Get metabolite fields

if numMets ~= numRxns % Distinguish using field lengths
    % Count number of elements in each field:
    nums = structfun(@numel, model);
    
    % Select fields with 'numRxns' elements:
    rxnFields = fieldNames(nums == numRxns);
else
    % Select all fields starting with 'rxn':
    rxnFields = fieldNames(strncmp('rxn', fieldNames, 3));
    
    % Select additional metabolite fields:
    rxnFields = [rxnFields; intersect(... % Additional reaction fields:
        {'lb', 'ub', 'rev', 'c', 'rules', 'grRules', 'confidenceScores', ...
        'ExchRxnBool', 'EXRxnBool', 'DMRxnBool', 'SinkRxnBool', 'SIntRxnBool'}, ...
        fieldNames)];
end

end
function metFields = GetMetFields(model)

%% Initialize

[numMets, numRxns] = size(model.S);

fieldNames = fieldnames(model);

%% Get metabolite fields

if numMets ~= numRxns % Distinguish using field lengths
    % Count number of elements in each field:
    nums = structfun(@numel, model);
    
    % Select fields with 'numMets' elements:
    metFields = fieldNames(nums == numMets);
else
    % Select all fields starting with 'met':
    metFields = fieldNames(strncmp('met', fieldNames, 3));
    
    % Select additional metabolite fields:
    metFields = [metFields; intersect(... % Additional metabolite fields:
        {'b'}, ...
        fieldNames)];
end

end
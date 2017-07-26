function subModel = ExtractSubModel(model, rxnInds)

metFields = GetMetFields(model); % Metabolite fields
rxnFields = GetRxnFields(model); % Reaction fields

metInvolved = any(model.S(:, rxnInds), 2); % Indices of involved metabolites

% Take part of stoichiometry matrix:
subModel.S = model.S(metInvolved, rxnInds);

% Add metabolite fields:
for iMetField = 1:numel(metFields)
    fieldName = metFields{iMetField}; % Current field
    
    vals = model.(fieldName); % Old field values
    
    subModel.(fieldName) = vals(metInvolved);
end

% Add reaction fields:
for iRxnField = 1:numel(rxnFields)
    fieldName = rxnFields{iRxnField}; % Current field
    
    vals = model.(fieldName); % Old field values
    
    subModel.(fieldName) = vals(rxnInds);
end

end
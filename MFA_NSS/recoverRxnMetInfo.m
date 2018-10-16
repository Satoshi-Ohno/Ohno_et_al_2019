function tmpModel = recoverRxnMetInfo(model)

tmpModel = model;

tmpFieldNames = fieldnames(tmpModel.rxnInfo);
for f = 1 : length(tmpFieldNames)
    tmpModel.(tmpFieldNames{f}) = tmpModel.rxnInfo.(tmpFieldNames{f});
end

tmpFieldNames = fieldnames(tmpModel.metInfo);
for f = 1 : length(tmpFieldNames)
    tmpModel.(tmpFieldNames{f}) = tmpModel.metInfo.(tmpFieldNames{f});
end

end
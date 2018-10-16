function field2var(inpustStruct)

tmpFieldNames = fieldnames(inpustStruct);
for i = 1 : length(tmpFieldNames)
    assignin('caller', tmpFieldNames{i}, inpustStruct.(tmpFieldNames{i}));
end

end
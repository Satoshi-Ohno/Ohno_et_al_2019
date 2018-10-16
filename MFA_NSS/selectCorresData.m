function selectedData = selectCorresData(data, idData, idSelect, dim)

idData = colVec(idData);
idSelect = colVec(idSelect);

tmpDiff = repmat(idData', length(idSelect), 1) - repmat(idSelect, 1, length(idData));
[~,idCorres] = min(abs(tmpDiff), [], 2);

switch dim
    case 1
        selectedData = data(idCorres, :);
    case 2
        selectedData = data(:, idCorres);
    case 3
        selectedData = data(:, :, idCorres);
end

end
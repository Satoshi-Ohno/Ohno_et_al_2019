function [maxmax, id] = mmax(data)

[maxVal, idRow] = max(data,[],1);
[maxmax, idCol] = max(maxVal);
id = [idRow(idCol), idCol];

end

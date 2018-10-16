function [minmin, id]= mmin(data)


[minVal, idRow] = min(data,[],1);
[minmin, idCol] = min(minVal);
id = [idRow(idCol), idCol];

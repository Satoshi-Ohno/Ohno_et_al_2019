%% Add constraints for lower and upper bounds of variables
function qpModel = modifyQPModelWithLbUb180116(model, expData, optionsMFA, qpModel, Cx, C0)
if isempty(qpModel.lb) 
    return
end

[nParam] = length(qpModel.q);

maxVal = max(abs([qpModel.lb, qpModel.ub]), [], 2);
isEqBound = qpModel.ub- qpModel.lb<=maxVal*0.01;
nEqBound = nnz(isEqBound);
nNonEqBound = nnz(~isEqBound);

tmpAeq = sparse(nEqBound, nParam);
tmpAeq(:,isEqBound) = tmpAeq(:,isEqBound)+ eye(nEqBound);

tmpA = sparse(nNonEqBound, nParam);
tmpA(:,~isEqBound) = tmpA(:,~isEqBound)+ eye(nNonEqBound);


qpModel.Aeq = [qpModel.Aeq; tmpAeq];
qpModel.beq = [qpModel.beq; (qpModel.lb(isEqBound) + qpModel.ub(isEqBound))/2];

qpModel.A = [qpModel.A; -tmpA; tmpA];
qpModel.b = [qpModel.b; -qpModel.lb(~isEqBound); qpModel.ub(~isEqBound)];


end

%% Add constrains for mass balance
function qpModel = modifyQPModelMassBalance(model, expData, optionsMFA, qpModel, optimType)
nConstr = size(optionsMFA.transformMat.constrMassBalance,1);
qpModel.Aeq = [qpModel.Aeq; optionsMFA.transformMat.constrMassBalance];
qpModel.beq = [qpModel.beq; zeros(nConstr,1)];

end

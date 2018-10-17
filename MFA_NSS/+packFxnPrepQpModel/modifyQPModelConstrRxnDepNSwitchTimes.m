%%  Add constrains for reaction dependent number of time intervals
function qpModel = modifyQPModelConstrRxnDepNKnots(model, expData, optionsMFA, qpModel, optimType)
tmpMatConstr=optionsMFA.convertMat.constrFluxNotForOptim;
nTmpConstr = size(tmpMatConstr,1);
tmpAeq1= tmpMatConstr;
tmpBeq1= zeros(nTmpConstr,1);

tmpMatConstr=optionsMFA.convertMat.constrConcRateNotForOptim;
nTmpConstr = size(tmpMatConstr,1);
tmpAeq2= tmpMatConstr;
tmpBeq2= zeros(nTmpConstr,1);

qpModel.Aeq = [qpModel.Aeq; tmpAeq1; tmpAeq2];
qpModel.beq = [qpModel.beq; tmpBeq1; tmpBeq2];

end

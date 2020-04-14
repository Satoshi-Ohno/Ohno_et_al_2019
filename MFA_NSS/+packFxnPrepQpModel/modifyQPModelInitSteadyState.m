%% Add constraints for steady state at initial time
function qpModel = modifyQPModelInitSteadyState(model, expData, optionsMFA, qpModel, optimType)

nSwitchTimes = optionsMFA.nSwitchTimes;
nNonPoolMets = optionsMFA.varSet.nNonPoolMets;
idNonPoolMets = optionsMFA.varSet.idNonPoolMets;

idParamInitConcRates = optionsMFA.idParamLocal.initConcRates;

isInitSteadyStateCompt = ismember(model.metInfo.compt(idNonPoolMets), optionsMFA.idComptInitSteadyState);

row = (1 : nnz(isInitSteadyStateCompt))';
col = idParamInitConcRates(isInitSteadyStateCompt);
tmpAeq = sparse(length(row), optionsMFA.idParamLocal.nParam);
tmpAeq(row, col) = speye(length(row), length(row));
tmpbeq = zeros(length(row),1);

qpModel.Aeq = [qpModel.Aeq; tmpAeq];
qpModel.beq = [qpModel.beq; tmpbeq];

end

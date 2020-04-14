%% Make transformation matrix (param -> conc)
function matParam2InitConcRate = ...
    makeMatParam2InitConcRate(model, expData, optionsMFA, fullSwitchTimes, inputFlag)

%% prepare variables
idParamLocal = optionsMFA.idParamLocal;
nParamLocal = idParamLocal.nParam;
% idParamMH = optionsMFA.idParamMH;
% nParamMH = idParamMH.nParam;
% nMets = optionsMFA.varSet.nMets;
% nRxns = optionsMFA.varSet.nRxns;
nNonPoolMets = optionsMFA.varSet.nNonPoolMets;
idNonPoolMets = optionsMFA.varSet.idNonPoolMets;
% nSwitchTimes = optionsMFA.varSet.nSwitchTimes;
% idEvalConcMets = optionsMFA.varSet.idEvalConcMets;
% nIndFluxes = optionsMFA.varSet.nIndFluxes;
% idNonIndFluxes = optionsMFA.varSet.idNonIndFluxes;
% nNonIndFluxes = optionsMFA.varSet.nNonIndFluxes;
% nNetRxns = optionsMFA.varSet.nNetRxns;
% nMHMetsConcRates = optionsMFA.varSet.nMHMetsConcRates;
% idMHMetsConcs = optionsMFA.varSet.idMHMetsConcs;
% idMHMetsConcRates = optionsMFA.varSet.idMHMetsConcRates;

%% make matrix
% inputFlag: 
% 11, isempty(fullSwitchTimes) == true & optionsMFA.isUseConcAsParam == true
% 12, isempty(fullSwitchTimes) == true & optionsMFA.isUseConcAsParam == false
% 21, isempty(fullSwitchTimes) == false & optionsMFA.isUseConcAsParam == true
% 22, isempty(fullSwitchTimes) == false & optionsMFA.isUseConcAsParam == false

switch inputFlag
    case {11,12}
        matParam2InitConcRate  = sparse(nNonPoolMets, nParamLocal);
        row = idNonPoolMets;
        col =  optionsMFA.idParamLocal.initConcRates;
        matParam2InitConcRate(row, col) = matParam2InitConcRate(row, col) + eye(nNonPoolMets);
    case {21,22}
        matParam2InitConcRate = optionsMFA.transformMat.param2InitConcRate;
end

end
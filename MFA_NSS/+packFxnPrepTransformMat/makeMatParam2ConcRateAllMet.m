%% Make transformation matrix (param -> concRateAllMet)
function matParam2ConcRateAllMet = ...
    makeMatParam2ConcRateAllMet(model, expData, optionsMFA, fullSwitchTimes, inputFlag)

%% prepare variables
% idParamLocal = optionsMFA.idParamLocal;
% nParamLocal = idParamLocal.nParam;
% idParamMH = optionsMFA.idParamMH;
% nParamMH = idParamMH.nParam;
nMets = optionsMFA.varSet.nMets;
nRxns = optionsMFA.varSet.nRxns;
% nNonPoolMets = optionsMFA.varSet.nNonPoolMets;
% idNonPoolMets = optionsMFA.varSet.idNonPoolMets;
nSwitchTimes = optionsMFA.varSet.nSwitchTimes;
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
    case 11
        matParam2ConcRateAllMet = [];
    case {12,21}
        matFlux2ConcRateAllMet = sparse(nMets*(nSwitchTimes+2), nRxns*(nSwitchTimes+2));
        tmpMat = full(model.S);
        for k = 1 : nSwitchTimes+2
            row = nMets*(k-1)+(1:nMets);
            col =  nRxns*(k-1)+(1:nRxns);
            matFlux2ConcRateAllMet(row, col) = matFlux2ConcRateAllMet(row, col)+tmpMat;
        end
        matParam2ConcRateAllMet = matFlux2ConcRateAllMet*optionsMFA.matParam2Flux;
    case 22
        matParam2ConcRateAllMet = optionsMFA.transformMat.param2ConcRateAllMet;
end

end
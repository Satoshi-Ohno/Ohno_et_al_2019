%% Make transformation matrix (paramMH -> concRate)
function matParamMH2ConcRate = ...
    makeMatParamMH2ConcRate(model, expData, optionsMFA, fullSwitchTimes, inputFlag)

%% prepare variables
% idParamLocal = optionsMFA.idParamLocal;
% nParamLocal = idParamLocal.nParam;
idParamMH = optionsMFA.idParamMH;
nParamMH = idParamMH.nParam;
% nMets = optionsMFA.varSet.nMets;
% nRxns = optionsMFA.varSet.nRxns;
% nNonPoolMets = optionsMFA.varSet.nNonPoolMets;
% idNonPoolMets = optionsMFA.varSet.idNonPoolMets;
nSwitchTimes = optionsMFA.varSet.nSwitchTimes;
% idEvalConcMets = optionsMFA.varSet.idEvalConcMets;
% nIndFluxes = optionsMFA.varSet.nIndFluxes;
% idNonIndFluxes = optionsMFA.varSet.idNonIndFluxes;
% nNonIndFluxes = optionsMFA.varSet.nNonIndFluxes;
% nNetRxns = optionsMFA.varSet.nNetRxns;
nMHMetsConcRates = optionsMFA.varSet.nMHMetsConcRates;
% idMHMetsConcs = optionsMFA.varSet.idMHMetsConcs;
% idMHMetsConcRates = optionsMFA.varSet.idMHMetsConcRates;

%% make matrix
% inputFlag: 
% 11, isempty(fullSwitchTimes) == true & optionsMFA.isUseConcAsParam == true
% 12, isempty(fullSwitchTimes) == true & optionsMFA.isUseConcAsParam == false
% 21, isempty(fullSwitchTimes) == false & optionsMFA.isUseConcAsParam == true
% 22, isempty(fullSwitchTimes) == false & optionsMFA.isUseConcAsParam == false

if ~optionsMFA.isUseQPInMH
    matParamMH2ConcRate = optionsMFA.matParam2ConcRate;
end

switch inputFlag
    case 11
        matParamMH2ConcRate = [];
    case 12
        matParamMH2ConcRate = sparse(nMHMetsConcRates*(nSwitchTimes+2), nParamMH);
        for k = 1 : nSwitchTimes+2
            row = nMHMetsConcRates*(k-1)+(1:nMHMetsConcRates);
            col =  idParamMH.concRates(nMHMetsConcRates*(k-1)+(1:nMHMetsConcRates));
            matParamMH2ConcRate(row, col) = matParamMH2ConcRate(row, col) + eye(nMHMetsConcRates);
        end
    case {21,22}
        matParamMH2ConcRate = optionsMFA.transformMat.paramMH2ConcRate;
end
    
end
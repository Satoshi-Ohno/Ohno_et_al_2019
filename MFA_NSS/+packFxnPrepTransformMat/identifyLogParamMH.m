%% Identify parameters in a log scale
function isLogParamMH= ...
    identifyLogParamMH(model, expData, optionsMFA, fullSwitchTimes, inputFlag)

%% prepare variables
idParamLocal = optionsMFA.idParamLocal;
nParamLocal = idParamLocal.nParam;
idParamMH = optionsMFA.idParamMH;
% nParamMH = idParamMH.nParam;
% nMets = optionsMFA.varSet.nMets;
% nRxns = optionsMFA.varSet.nRxns;
% nNonPoolMets = optionsMFA.varSet.nNonPoolMets;
% idNonPoolMets = optionsMFA.varSet.idNonPoolMets;
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
%         keyboard;
        matParamMH2ParamLocal = ...
            packFxnPrepTransformMat.makeMatParamMH2ParamLocal(model, expData, optionsMFA, fullSwitchTimes, inputFlag);        
        matParamLocal2ParamMH = matParamMH2ParamLocal';

        isLogParamLocal = false(nParamLocal,1);
        idLog = [idParamLocal.concs, idParamLocal.fluxes, idParamLocal.switchTimes];
        isLogParamLocal(idLog) = true;
        
        isLogParamMH = logical(matParamLocal2ParamMH*isLogParamLocal);

    case {21,22}
        isLogParamMH = optionsMFA.transformMat.isLogParamMH;
end
end
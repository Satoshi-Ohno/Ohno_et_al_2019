%% Make transformation matrix (paramMH -> ind)
function matParamMH2Ind = ...
    makeMatParamMH2Ind(model, expData, optionsMFA, fullSwitchTimes, inputFlag)

%% prepare variables
% idParamLocal = optionsMFA.idParamLocal;
% nParamLocal = idParamLocal.nParam;
idParamMH = optionsMFA.idParamMH;
nParamMH = idParamMH.nParam;
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
        if isfield(optionsMFA, 'isIndParams') && isfield(optionsMFA.isIndParams, 'MH')
            idIndParamsMH = find(optionsMFA.isIndParams.MH);
            matParamMH2Ind = sparse(length(idIndParamsMH), nParamMH);
            matParamMH2Ind(:, idIndParamsMH) = speye(length(idIndParamsMH));
        else
            matParamMH2Ind = [];
        end
    case {21,22}
        matParamMH2Ind = optionsMFA.transformMat.paramMH2Ind;
end
    
end
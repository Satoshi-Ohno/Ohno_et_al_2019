%% Make transformation matrix (paramMH -> paramLocal)
function matParamMH2ParamLocal = ...
    makeMatParamMH2ParamLocal(model, expData, optionsMFA, fullSwitchTimes, inputFlag)

%% prepare variables
idParamLocal = optionsMFA.idParamLocal;
nParamLocal = idParamLocal.nParam;
idParamMH = optionsMFA.idParamMH;
nParamMH = idParamMH.nParam;
% nMets = optionsMFA.varSet.nMets;
% nRxns = optionsMFA.varSet.nRxns;
nNonPoolMets = optionsMFA.varSet.nNonPoolMets;
% idNonPoolMets = optionsMFA.varSet.idNonPoolMets;
nSwitchTimes = optionsMFA.varSet.nSwitchTimes;
% idEvalConcMets = optionsMFA.varSet.idEvalConcMets;
% nIndFluxes = optionsMFA.varSet.nIndFluxes;
% idNonIndFluxes = optionsMFA.varSet.idNonIndFluxes;
% nNonIndFluxes = optionsMFA.varSet.nNonIndFluxes;
% nNetRxns = optionsMFA.varSet.nNetRxns;
% nMHMetsConcRates = optionsMFA.varSet.nMHMetsConcRates;
idMHMetsConcs = optionsMFA.varSet.idMHMetsConcs;
idMHMetsConcRates = optionsMFA.varSet.idMHMetsConcRates;

%% make matrix
% inputFlag: 
% 11, isempty(fullSwitchTimes) == true & optionsMFA.isUseConcAsParam == true
% 12, isempty(fullSwitchTimes) == true & optionsMFA.isUseConcAsParam == false
% 21, isempty(fullSwitchTimes) == false & optionsMFA.isUseConcAsParam == true
% 22, isempty(fullSwitchTimes) == false & optionsMFA.isUseConcAsParam == false

if ~optionsMFA.isUseQPInMH
    matParamMH2ParamLocal = speye(nParamLocal, nParamMH);
    return
end

switch inputFlag
    case 11
        matParamMH2ParamLocal = sparse(nParamLocal, nParamMH);
        tmpFieldNames = fieldnames(idParamMH);
        for f = 1 : length(tmpFieldNames)
            switch tmpFieldNames{f}
                case {'fluxes', 'switchTimes'}
                    row = idParamLocal.(tmpFieldNames{f});
                    col = idParamMH.(tmpFieldNames{f});
                case {'initConcRates'}
                    row = idParamLocal.(tmpFieldNames{f})(idMHMetsConcRates);
                    col = idParamMH.(tmpFieldNames{f});
                case {'concs'}
                    isConcMHMets = false(nNonPoolMets, nSwitchTimes+2);
                    isConcMHMets(idMHMetsConcs,:) = true;
                    isConcMHMets = reshape(isConcMHMets, nNonPoolMets*(nSwitchTimes+2), 1);
                    row = idParamLocal.concs(isConcMHMets);
                    col = idParamMH.(tmpFieldNames{f});
                otherwise
                    row = [];
                    col = [];
            end
            matParamMH2ParamLocal(row, col) = matParamMH2ParamLocal(row, col) + eye(length(row));
        end
    case 12
        matParamMH2ParamLocal = sparse(nParamLocal, nParamMH);
        tmpFieldNames = fieldnames(idParamMH);
        for f = 1 : length(tmpFieldNames)
            switch tmpFieldNames{f}
                case {'fluxes', 'switchTimes'}
                    row = idParamLocal.(tmpFieldNames{f});
                    col = idParamMH.(tmpFieldNames{f});
                case {'initConcs'}
                    row = idParamLocal.(tmpFieldNames{f})(idMHMetsConcs);
                    col = idParamMH.(tmpFieldNames{f});
                case {'concRates'}
                    isConcRateOuterOptimMets = false(nNonPoolMets, nSwitchTimes+2);
                    isConcRateOuterOptimMets(idMHMetsConcRates,:) = true;
                    isConcRateOuterOptimMets = reshape(isConcRateOuterOptimMets, nNonPoolMets*(nSwitchTimes+2), 1);
                    row = idParamLocal.concRates(isConcRateOuterOptimMets);
                    col = idParamMH.(tmpFieldNames{f});
                otherwise
                    row = [];
                    col = [];
            end
            matParamMH2ParamLocal(row, col) = matParamMH2ParamLocal(row, col) + eye(length(row));
        end
    case {21,22}
        matParamMH2ParamLocal = optionsMFA.transformMat.paramMH2ParamLocal;
end
    
end
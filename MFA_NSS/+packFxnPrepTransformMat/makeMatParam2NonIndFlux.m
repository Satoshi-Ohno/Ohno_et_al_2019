%% Make transformation matrix (param -> netFlux)
function matParam2NonIndFlux = ...
    makeMatParam2NonIndFlux(model, expData, optionsMFA, fullSwitchTimes, inputFlag)

%% prepare variables
% idParamLocal = optionsMFA.idParamLocal;
% nParamLocal = idParamLocal.nParam;
% idParamMH = optionsMFA.idParamMH;
% nParamMH = idParamMH.nParam;
% nMets = optionsMFA.varSet.nMets;
nRxns = optionsMFA.varSet.nRxns;
% nNonPoolMets = optionsMFA.varSet.nNonPoolMets;
% idNonPoolMets = optionsMFA.varSet.idNonPoolMets;
nSwitchTimes = optionsMFA.varSet.nSwitchTimes;
% idEvalConcMets = optionsMFA.varSet.idEvalConcMets;
% nIndFluxes = optionsMFA.varSet.nIndFluxes;
idNonIndFluxes = optionsMFA.varSet.idNonIndFluxes;
nNonIndFluxes = optionsMFA.varSet.nNonIndFluxes;
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
        matParam2NonIndFlux = [];
    case {12,21}
        matFlux2NonIndFlux = sparse(nNonIndFluxes*(nSwitchTimes+2), nRxns*(nSwitchTimes+2));
        for k = 1 : nSwitchTimes+2
            row = nNonIndFluxes*(k-1)+(1:nNonIndFluxes);
            col =  nRxns*(k-1)+idNonIndFluxes;
            matFlux2NonIndFlux(row, col) = matFlux2NonIndFlux(row, col) + eye(nNonIndFluxes);
        end        
        matParam2NonIndFlux = matFlux2NonIndFlux*optionsMFA.matParam2Flux;
    case 22
        matParam2NonIndFlux= optionsMFA.transformMat.param2NonIndFlux;
end

end
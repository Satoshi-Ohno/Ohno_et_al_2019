%% Make transformation matrix (param -> conc)
function matParam2Conc = ...
    makeMatParam2Conc(model, expData, optionsMFA, fullSwitchTimes, inputFlag)

%% prepare variables
idParamLocal = optionsMFA.idParamLocal;
nParamLocal = idParamLocal.nParam;
% idParamMH = optionsMFA.idParamMH;
% nParamMH = idParamMH.nParam;
% nMets = optionsMFA.varSet.nMets;
% nRxns = optionsMFA.varSet.nRxns;
nNonPoolMets = optionsMFA.varSet.nNonPoolMets;
idNonPoolMets = optionsMFA.varSet.idNonPoolMets;
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
        matParam2Conc = sparse(nNonPoolMets*(nSwitchTimes+2), nParamLocal);
        for k = 1 : nSwitchTimes+2
            row = nNonPoolMets*(k-1)+(1:nNonPoolMets);
            col =  idParamLocal.concs(nNonPoolMets*(k-1)+(1:nNonPoolMets));
            matParam2Conc(row, col) = matParam2Conc(row, col) + eye(nNonPoolMets);
        end
    case 12
        matParam2Conc = [];
    case 21
        matParam2Conc = optionsMFA.transformMat.param2Conc;
    case 22
        matParam2ConcConst = optionsMFA.matParam2ConcConst;
        matParam2ConcRate = optionsMFA.transformMat.param2ConcRate;
        matConcConstRate2Conc = ...
            packFxnPrepTransformMat.makeMatConcConstRate2Conc(model, expData, optionsMFA, fullSwitchTimes, fullSwitchTimes, idNonPoolMets);
        matParam2Conc = matConcConstRate2Conc * [matParam2ConcConst;matParam2ConcRate];
end

end
%% Make transformation matrix (param -> concExpTime)
function matParam2ConcExpTime = ...
    makeMatParam2ConcExpTime(model, expData, optionsMFA, fullSwitchTimes, inputFlag)

import packFxnPrepTransformMat.makeMatConcConstRate2Conc

%% prepare variables
% idParamLocal = optionsMFA.idParamLocal;
% nParamLocal = idParamLocal.nParam;
% idParamMH = optionsMFA.idParamMH;
% nParamMH = idParamMH.nParam;
% nMets = optionsMFA.varSet.nMets;
% nRxns = optionsMFA.varSet.nRxns;
% nNonPoolMets = optionsMFA.varSet.nNonPoolMets;
% idNonPoolMets = optionsMFA.varSet.idNonPoolMets;
% nSwitchTimes = optionsMFA.varSet.nSwitchTimes;
idEvalConcMets = optionsMFA.varSet.idEvalConcMets;
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
    case {11, 12}
        matParam2ConcExpTime = [];
    case {21,22}
        matParam2ConcConst = optionsMFA.matParam2ConcConst;
        matParam2ConcRate = optionsMFA.matParam2ConcRate;
        matConcConstRate2ConcExpTime = ...
            makeMatConcConstRate2Conc(model, expData, optionsMFA, fullSwitchTimes, expData.time, idEvalConcMets);
        matParam2ConcExpTime = matConcConstRate2ConcExpTime * [matParam2ConcConst;matParam2ConcRate];
end

end
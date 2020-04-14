%% Make matrix to constrain concRate which is not used fof optimization
function matConstrFluxNotForOptim = ...
    makeMatConstrFluxNotForOpim(model, expData, optionsMFA, fullSwitchTimes, inputFlag)

import packFxnPrepTransformMat.makeVar2NotForOptim
import packFxnPrepTransformMat.makeCalcVar2NotForOptim

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
    case {11, 12}
        matConstrFluxNotForOptim = [];
    case {21,22}
        varType = 'flux';
        matFlux2NotForOptim = ...
            makeVar2NotForOptim(model, expData, optionsMFA, fullSwitchTimes, varType);
        matCalcFlux2NotForOptim = ...
            makeCalcVar2NotForOptim(model, expData, optionsMFA, fullSwitchTimes, varType);
        matConstrFluxNotForOptim = ...
            (matCalcFlux2NotForOptim - matFlux2NotForOptim) * optionsMFA.matParam2Flux;        
end

end



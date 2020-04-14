%% Make transformation matrix to obtain variales from parameter vector and switch time
function transformMat = prepTransformMat(model, expData, optionsMFA, fullSwitchTimes)

if nargin <=3
    fullSwitchTimes = [];
end

import packFxnPrepTransformMat.*

%% input flag
if isempty(fullSwitchTimes)
    if optionsMFA.isUseConcAsParam
        inputFlag = 11;
    else
        inputFlag = 12;
    end
else
    if optionsMFA.isUseConcAsParam
        inputFlag = 21;
    else
        inputFlag = 22;
    end
end

%% Make transformation matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% param -> concRate
matParam2ConcRate = ...
    makeMatParam2ConcRate(model, expData, optionsMFA, fullSwitchTimes, inputFlag);
optionsMFA.matParam2ConcRate = matParam2ConcRate;

%% param -> concConst
matParam2ConcConst =  ...
    makeMatParam2ConcConst(model, expData, optionsMFA, fullSwitchTimes, inputFlag);
optionsMFA.matParam2ConcConst = matParam2ConcConst;

%% param -> initConc
matParam2InitConc = ...
    makeMatParam2InitConc(model, expData, optionsMFA, fullSwitchTimes, inputFlag);

%% param -> conc
matParam2Conc = ...
    makeMatParam2Conc(model, expData, optionsMFA, fullSwitchTimes, inputFlag);

%% param -> initConcRate
matParam2InitConcRate = ...
    makeMatParam2InitConcRate(model, expData, optionsMFA, fullSwitchTimes, inputFlag);

%% param -> indFlux
matParam2IndFlux = ...
    makeMatParam2IndFlux(model, expData, optionsMFA, fullSwitchTimes, inputFlag);
optionsMFA.matParam2IndFlux = matParam2IndFlux;

%% param -> flux
matParam2Flux = ...
    makeMatParam2Flux(model, expData, optionsMFA, fullSwitchTimes, inputFlag);
optionsMFA.matParam2Flux = matParam2Flux;

%% param -> nonIndFlux
matParam2NonIndFlux = ...
    makeMatParam2NonIndFlux(model, expData, optionsMFA, fullSwitchTimes, inputFlag);
optionsMFA.matParam2NonIndFlux = matParam2NonIndFlux;

%% param -> netFlux
matParam2NetFlux = ...
    makeMatParam2NetFlux(model, expData, optionsMFA, fullSwitchTimes, inputFlag);

%% param -> concRateAllMet
matParam2ConcRateAllMet = ...
    makeMatParam2ConcRateAllMet(model, expData, optionsMFA, fullSwitchTimes, inputFlag);

%% paramMH -> paramLocal
matParamMH2ParamLocal = ...
    makeMatParamMH2ParamLocal(model, expData, optionsMFA, fullSwitchTimes, inputFlag);

%% paramLocal -> paramMH
matParamLocal2ParamMH = matParamMH2ParamLocal';

%% Identify parameters in a log scale
isLogParamMH= ...
    identifyLogParamMH(model, expData, optionsMFA, fullSwitchTimes, inputFlag);

%% paramMH -> ind
matParamMH2Ind = ...
    makeMatParamMH2Ind(model, expData, optionsMFA, fullSwitchTimes, inputFlag);

%% paramMH -> concRate
% This may not be necessary
matParamMH2ConcRate = ...
    makeMatParamMH2ConcRate(model, expData, optionsMFA, fullSwitchTimes, inputFlag);

%% paramMH -> indFlux
% This may not be necessary
matParamMH2IndFlux = ...
    makeMatParamMH2IndFlux(model, expData, optionsMFA, fullSwitchTimes, inputFlag);

%% constraint for mass balance
matConstrMassBalance = ...
    makeMatConstrMassBalance(model, expData, optionsMFA, fullSwitchTimes, inputFlag);

%% param -> concQP
matParam2ConcQP = ...
    makeMatParam2ConcQP(model, expData, optionsMFA, fullSwitchTimes, inputFlag);

%% param -> concExpTime (only idEvalConcMets)
matParam2ConcExpTime = ...
    makeMatParam2ConcExpTime(model, expData, optionsMFA, fullSwitchTimes, inputFlag);

%% param -> pCoefConc
matParam2pCoefConc = ...
    makeMatParam2pCoefConc(model, expData, optionsMFA, fullSwitchTimes, inputFlag);

%% param -> pCoefFlux
matParam2pCoefFlux = ...
    makeMatParam2pCoefFlux(model, expData, optionsMFA, fullSwitchTimes, inputFlag);

%% Constraint for reaction dependent number of time intervvals
matConstrConcRateNotForOptim = ...
    makeMatConstrConcRateNotForOpim(model, expData, optionsMFA, fullSwitchTimes, inputFlag);
matConstrFluxNotForOptim = ...
    makeMatConstrFluxNotForOpim(model, expData, optionsMFA, fullSwitchTimes, inputFlag);


%% Make structure of transformation matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%
transformMat.param2InitConc = matParam2InitConc;
transformMat.param2Conc = matParam2Conc;
transformMat.param2InitConcRate = matParam2InitConcRate;
transformMat.param2ConcRate = matParam2ConcRate;
transformMat.param2ConcRateAllMet = matParam2ConcRateAllMet;
transformMat.param2IndFlux = matParam2IndFlux;
transformMat.param2Flux = matParam2Flux;
transformMat.param2NonIndFlux = matParam2NonIndFlux;
transformMat.param2NetFlux = matParam2NetFlux;
transformMat.paramMH2ParamLocal = matParamMH2ParamLocal;
transformMat.paramLocal2ParamMH = matParamLocal2ParamMH;
transformMat.paramMH2Ind = matParamMH2Ind;
transformMat.paramMH2ConcRate = matParamMH2ConcRate;
transformMat.paramMH2IndFlux = matParamMH2IndFlux;
transformMat.constrMassBalance = matConstrMassBalance;

tmpFieldNames = fieldnames(transformMat);
for f = 1 : length(tmpFieldNames)
    tmpMat = transformMat.(tmpFieldNames{f});
    tmpMat(abs(tmpMat)<=10^-9) = 0;  % remove error
    transformMat.(tmpFieldNames{f})= tmpMat;
end

transformMat.param2ConcConst = matParam2ConcConst;
transformMat.param2ConcQP = matParam2ConcQP;
transformMat.param2ConcExpTime = matParam2ConcExpTime;
transformMat.param2pCoefConc = matParam2pCoefConc;
transformMat.param2pCoefFlux = matParam2pCoefFlux;
transformMat.constrConcRateNotForOptim = matConstrConcRateNotForOptim;
transformMat.constrFluxNotForOptim = matConstrFluxNotForOptim;

tmpFieldNames = fieldnames(transformMat);
for f = 1 : length(tmpFieldNames)
    tmpMat = transformMat.(tmpFieldNames{f});
    tmpMat(abs(tmpMat)<=10^-9) = 0;  % remove error
    transformMat.(tmpFieldNames{f})= tmpMat;
end

transformMat.isLogParamMH = isLogParamMH;

end

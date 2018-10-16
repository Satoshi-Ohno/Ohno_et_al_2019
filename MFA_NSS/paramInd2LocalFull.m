function [paramLocal, convertMat] = paramInd2LocalFull180104(model, expData, optionsMFA, paramInd, optimType)

%% Make QP model

isInputIndVars = true;
isInputFullIndVarsConstr= true;
fxnPrepQpModel = str2func(optionsMFA.fxnName.prepQpModel);
if optionsMFA.isUseQPInMH
    [qpModel, convertMat] = fxnPrepQpModel(model, expData, optionsMFA,  ...
        optimType, paramInd,isInputIndVars, isInputFullIndVarsConstr);
else
    switch optimType
        case {'local'}
            [qpModel, convertMat] = fxnPrepQpModel(model, expData, optionsMFA,  ...
                optimType, paramInd,isInputIndVars, isInputFullIndVarsConstr);
        case {'metaheuristic', 'init', 'QP'}
            tmpOptimType = 'local';
            tmpParamInd = paramInd;
            isIndParam = optionsMFA.isIndParams.MH;
            tmpParamInd(optionsMFA.convertMat.isLogParamMH(isIndParam)) = 10.^paramInd(optionsMFA.convertMat.isLogParamMH(isIndParam));
            [qpModel, convertMat] = fxnPrepQpModel(model, expData, optionsMFA,  ...
                tmpOptimType, tmpParamInd,isInputIndVars, isInputFullIndVarsConstr);
    end
end

%% Calculate paramLocal using QP model
if optionsMFA.isUseQPInMH
    switch optimType
        case {'local'}
            paramLocal = qpModel.L*paramInd + qpModel.m;
        case {'metaheuristic', 'init', 'QP'}
                isMinVar = true;
                fxnModifyQPModelObjFun = str2func(optionsMFA.fxnName.modifyQPModelObjFun);
                qpModel = fxnModifyQPModelObjFun(model, expData, optionsMFA, qpModel, isMinVar);
                [qpSol, qpScore, exitFlag] = solveQP(qpModel);
                if isempty(qpSol)
                    paramLocal = [];
                    return
                end
                paramLocal = qpModel.L*qpSol + qpModel.m;
    end
else
    switch optimType
        case {'local'}
            paramLocal = qpModel.L*paramInd + qpModel.m;
        case {'metaheuristic', 'init', 'QP'}
            isIndParam = optionsMFA.isIndParams.MH;
            paramInd(optionsMFA.convertMat.isLogParamMH(isIndParam)) = 10.^paramInd(optionsMFA.convertMat.isLogParamMH(isIndParam));
            paramLocal = qpModel.L*paramInd + qpModel.m;
    end
end

end
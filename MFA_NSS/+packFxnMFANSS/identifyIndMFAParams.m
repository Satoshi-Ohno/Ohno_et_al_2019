%% Identify independent parameters
function [isIndMFAParams, isRedundantConstrAeq] = ...
    identifyIndMFAParams(model, expData, optionsMFA, optimType)

field2var(optionsMFA.varSet)
idParamLocal = optionsMFA.idParamLocal;
nParamLocal = idParamLocal.nParam;
lbInit = optionsMFA.lbInit;
ubInit = optionsMFA.ubInit;
optimTypeInit = 'init';
idParamMH = optionsMFA.idParamMH;
nParamMH = idParamMH.nParam;

%% temporary parameter
tmpParamMH = rand(nParamMH, 1)*3;

%% Initial switch time
log10Range = log10(ubInit.switchTimes)-log10(lbInit.switchTimes);
while true
    switchTimes = 10.^(log10Range.*rand(nSwitchTimes, 1)+log10(lbInit.switchTimes));
    switchTimes = sort(switchTimes);
    fullSwitchTimes = [0, switchTimes', expData.time(end)];
    if all(diff(fullSwitchTimes) >= optionsMFA.minKnotTimeDiff)
        break
    end
end

tmpParamMH(idParamMH.switchTimes) = log10(switchTimes);

%% Convert paramMH to independent parameters
switch optimType
    case {'QP'}
        tmpParamMH = optionsMFA.transformMat.paramMH2Ind * tmpParamMH;
end

%% Make QP model
switch optimType
    case {'QP'}
        isInputIndVars = true;
    otherwise
        isInputIndVars = false;
end
isInputFullIndVarsConstr= false;
[qpModel, ~, qpModelConstr] = prepQpModel(model, expData, optionsMFA,  ...
    optimType, tmpParamMH,isInputIndVars, isInputFullIndVarsConstr);

%% Modify objective function of QP
isMinVar = true;
qpModel = modifyQPModelObjFun(model, expData, optionsMFA, qpModel, isMinVar);

%% Solve QP to fit metabolite concentrations
[qpSol, qpScore,exitFlag] = solveQP(qpModel);
isIndMFAParams = false(nParamLocal,1);
isIndMFAParams(qpModel.idIndVars) = true;
switch optimType
    case {'init', 'local', 'QP'}
    case {'metaheuristic'}
        isIndMFAParams = optionsMFA.transformMat.paramLocal2ParamMH * isIndMFAParams==true;
end
isRedundantConstrAeq = false(size(qpModelConstr.Aeq,1), 1);
isRedundantConstrAeq(qpModel.idRedundantConstr) = true;

end

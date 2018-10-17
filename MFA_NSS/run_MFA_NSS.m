function [sol, model, expData, optionsMFA] =  run_MFA_NSS(xlsModelFileName, expData, optionsMFA)
%run_MFA_NSS        run metabolic flux analysis under non-steady-state
% INPUT
% xlsModelFileName        Excel file name for model
% expData         Measured data whch is used to calculate wRSS
% optionsMFA          options of metaboic flux analysis under non-steady-state
% 
% OUTPUT
% sol         solution of metabolic flux analysis under non-steady-state
% model       metabolic network model
% expData       expData wher some fields are added
% optionsMFA        optionsMFA where some fields are added

import packFxnMIMFA.*

%% Load default options
defaultOptionsMFA = loadDefaultOptionsMFA();
filedNamesInputOptionsMFA = fieldnames(optionsMFA);
for f = 1 : length(filedNamesInputOptionsMFA)
    defaultOptionsMFA.(filedNamesInputOptionsMFA{f}) = optionsMFA.(filedNamesInputOptionsMFA{f});
end
optionsMFA = defaultOptionsMFA;

if optionsMFA.isRxnDepNSwitchTimes
    optionsMFA.nSwitchTimes = max(optionsMFA.nSwitchTimesNetRxns);
    optionsMFA.patternNSwitchTimes = unique(optionsMFA.nSwitchTimesNetRxns);
else
    optionsMFA.patternNSwitchTimes = optionsMFA.nSwitchTimes;
end

%% load model defined in Excel file
model.name = strrep(xlsModelFileName, '.xls', '');
model = loadXlsNetwork(model, optionsMFA, xlsModelFileName);

%% Add information on model
% Stoichiometic matrix
model = makeSMatrix(model, optionsMFA);

% EMU network
model = createEMUs(model, optionsMFA);

% EMU matrix
model = createEmuMatrix(model, optionsMFA);

% Matrix for calcuation of net flux from forwad and backward flux
model = makeNetFluxInfo(model, optionsMFA);


%% Identification of independent flux
A = full(model.S(~model.metInfo.isPoolMets, :));
idInputIndVars = [];
isInputFullIndVarsConstr = false;
fxnSelectIndVars = str2func(optionsMFA.fxnName.selectIndVars);
[idIndVars, idRedandantConstr, Cx, C0, idInvalidIndFluxes] = ...
    fxnSelectIndVars(A, isInputFullIndVarsConstr, idInputIndVars);
model.idIndFluxes = idIndVars;
model.invS = [C0,Cx];

disp('Following flux cannot be independent. Removed from input independent flux.')
for i = 1 :length(idInvalidIndFluxes)
    disp(['# ' num2str(idInvalidIndFluxes(i)) ': ' model.rxnNames{idInvalidIndFluxes(i)} ''])
end

%% Prepare optimization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varSet = prepVarsMFA(model, expData, optionsMFA);
optionsMFA.varSet = varSet;

optionsMFA = prepNSwitchTimesRxnsMets(model, expData, optionsMFA);

[optionsMFA.lbInit, optionsMFA.ubInit] = ...
    prepParamBound(model, expData, optionsMFA, 'init', 1/optionsMFA.paramBoundExpandRatio.^2);
[optionsMFA.lbMH, optionsMFA.ubMH] = ...
    prepParamBound(model, expData, optionsMFA, 'metaheuristic', 1/optionsMFA.paramBoundExpandRatio);
[optionsMFA.lbLocal, optionsMFA.ubLocal] = ...
    prepParamBound(model, expData, optionsMFA, 'local');

optionsMFA.idParamMH= ...
    prepIdParam(model, expData, optionsMFA, 'metaheuristic');
optionsMFA.idParamLocal = ...
    prepIdParam(model, expData, optionsMFA, 'local');

fxnPrepConvertMat = str2func(optionsMFA.fxnName.prepConvertMat);
optionsMFA.convertMat = fxnPrepConvertMat(model, expData, optionsMFA);

expData.vecExpConcStruct = makeVecExpConcs(model, expData, optionsMFA);

%% Identify independent parameters
[optionsMFA.isIndParams.init, optionsMFA.isRedundantConstrAeq.init] = ...
    identifyIndMFAParams(model, expData, optionsMFA, 'init');
[optionsMFA.isIndParams.MH, optionsMFA.isRedundantConstrAeq.MH] = ...
    identifyIndMFAParams(model, expData, optionsMFA, 'metaheuristic');
[optionsMFA.isIndParams.local, optionsMFA.isRedundantConstrAeq.local] = ...
    identifyIndMFAParams(model, expData, optionsMFA, 'local');

optionsMFA.convertMat = fxnPrepConvertMat(model, expData, optionsMFA);
if optionsMFA.isUseQPInMH 
    [optionsMFA.isIndParams.QP, optionsMFA.isRedundantConstrAeq.QP] = ...
        identifyIndMFAParams(model, expData, optionsMFA, 'QP');
end

%% Options for ODE 
for i = 1 : length(optionsMFA.optionsODE.activeOptionNames);
    tmpOptionName = optionsMFA.optionsODE.activeOptionNames{i};
    if i == 1
        optionsMFA.optionsODE.final  = odeset(tmpOptionName, optionsMFA.optionsODE.(tmpOptionName));
    else
        optionsMFA.optionsODE.final= odeset(optionsMFA.optionsODE.final, tmpOptionName, optionsMFA.optionsODE.(tmpOptionName));
    end
end

disp('MI-MFA start')
if optionsMFA.isRxnDepNSwitchTimes
    disp(['reaction-specific # of time intervals'])
else
    disp(['# of time intervals: ' num2str(optionsMFA.nSwitchTimes+1) ...
        ' (' num2str(optionsMFA.nSwitchTimes) ' switch times)'])
end


%% Initial parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fxnSolveMDVDynamics = str2func(optionsMFA.fxnName.solveMDVDynamics);
fxnCreateInitParams = str2func(optionsMFA.fxnName.createInitParams);
[initParams, initScores, nExpData] =  fxnCreateInitParams(...
    model, expData,optionsMFA);

expData.nExpData = nExpData;
optionsMFA.optionsOptimMH.initParams = initParams;
optionsMFA.optionsOptimMH.initScores = initScores;

optimType =  'metaheuristic';
objfun =@(params) fxnSolveMDVDynamics(...
    params, model, expData, optionsMFA, optimType);
[~,~,sol] = objfun(initParams(:,1));

disp(['Best Score by initial sampling before optimization:' num2str(min(initScores))])




%% Metaheuristic optimization %%%%%%%%%%%%%%%%%%%%%%%%%%%%
optimType =  'metaheuristic';
objfun =@(params) fxnSolveMDVDynamics(...
    params, model, expData, optionsMFA, optimType);

disp(' ')
disp('Metaheuristic optimization by CMA-ES')
solMH = performOptimMH(model, expData, optionsMFA, objfun, initParams);
disp(['Best Score by metaheuristic optimization:' num2str(solMH.score)])

%% Local optimization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
optimType =  'local';
objfun =@(params) fxnSolveMDVDynamics(...
    params, model, expData, optionsMFA, optimType);

initParamLocal = solMH.paramVector;
optionsMFA.optionsOptimLocal.initParams = initParamLocal;

disp(' ')
disp('Local optimization.')
solLocal = performOptimLocal(model, expData, optionsMFA, objfun, initParamLocal);
if ~isempty(solLocal)
    disp(['Best Score by local optimization:' num2str(solLocal.score)])
end

%% Final solution
if isempty(solLocal)
    sol = solMH;
    sol.optimType = 'metaheuristic';
else
    if solLocal.score <= solMH.score*(1+sign(solMH.score)*10^-6)
        sol = solLocal;
        sol.optimType = 'local';
    else
        sol = solMH;
        sol.optimType = 'metaheuristic';
    end
end

end

%% Metaheuristic optimization
function  solMH = performOptimMH(model, expData, optionsMFA, objfun, initParams)

optionsOptimMH = optionsMFA.optionsOptimMH;
optimType =  'metaheuristic';
isLogTransform = true;

%% Set parameter bound
isIndParam = true;
lbMH = optionsMFA.lbMH;
ubMH = optionsMFA.ubMH;
fxnParam2vector = str2func(optionsMFA.fxnName.param2vector);
optionsOptimMH.lb = fxnParam2vector(lbMH, optionsMFA, optimType, isLogTransform, isIndParam);
optionsOptimMH.ub = fxnParam2vector(ubMH, optionsMFA, optimType, isLogTransform, isIndParam);


%% Setting options of CMA-ES
tmpSigma = sqrt(var(initParams,[],2));

objfunStr = optionsMFA.fxnName.solveMDVDynamics;
fxnSolveMDVDynamics = str2func(objfunStr);
objfun =@(params) fxnSolveMDVDynamics(...
    params, model, expData, optionsMFA, optimType);
sigmaCMAES = tmpSigma;
varaginCMAES{1} = model;
varaginCMAES{2} = expData;
varaginCMAES{3} = optionsMFA;
varaginCMAES{4} = optimType;
optionsCMAES = cmaes('defaults');
optionsCMAES.MaxFunEvals = optionsOptimMH.nMaxEval;
optionsCMAES.LBounds = optionsOptimMH.lb;
optionsCMAES.UBounds = optionsOptimMH.ub;
optionsCMAES.PopSize = optionsOptimMH.popSize;
optionsCMAES.DispModulo = ceil(optionsOptimMH.nMaxEval*0.01/(optionsOptimMH.popSize/2)); % nMaxEval‚Ì1%‚¸‚ÂB
if ~ isempty(optionsOptimMH.CMA.active)
    optionsCMAES.CMA.active = optionsOptimMH.CMA.active;
end

%% Run metaheuristic optimization
[solMHRaw.bestParam, solMHRaw.bestScore, ...
    solMHRaw.couteval, solMHRaw.stopflag,solMHRaw.out,solMHRaw.bestever] = ...
    cmaes(objfunStr, initParams, sigmaCMAES, optionsCMAES, varaginCMAES{:});

%% Final solution of metaheuristic optimization
[~,~,solMH] = objfun(solMHRaw.bestever.x);

end

%% Local optimization
function solLocal = performOptimLocal(model, expData, optionsMFA, objfun, initParamLocal)
optionsOptimLocal = optionsMFA.optionsOptimLocal;
optimType = 'local';
isLogTransform = false;
isIndParam = true;

%% Set parameter bounds
lbLocal = optionsMFA.lbLocal;
ubLocal = optionsMFA.ubLocal;

fxnParam2vector = str2func(optionsMFA.fxnName.param2vector);
optionsOptimLocal.lb = fxnParam2vector(lbLocal, optionsMFA, optimType, isLogTransform, isIndParam);
optionsOptimLocal.ub = fxnParam2vector(ubLocal, optionsMFA, optimType, isLogTransform, isIndParam);

%% Initial parameters
initParamLocal = initParamLocal(optionsMFA.isIndParams.local);

%% Confirm parameters in the bounds
if any(initParamLocal < optionsOptimLocal.lb)
    find(initParamLocal < optionsOptimLocal.lb)'
    error('Some of initial parameters are smaller than lower bound in local optimization.')
end
if any(initParamLocal > optionsOptimLocal.ub)
    find(initParamLocal> optionsOptimLocal.ub)'
    error('Some of initial parameters are larger than upper bound in local optimization.')
end

%% Define function to calculate residuals
objfunRes = @(params) calcRes(objfun, params);

%% Set constrains
[qpModel] = makeQpModelFmincon(model, expData, optionsMFA, initParamLocal);
A = qpModel.A;
b = qpModel.b;
Aeq = qpModel.Aeq;
beq = qpModel.beq;
nonlcon = @(paramLocal) makeNonLinConstrFmincon(model, expData, optionsMFA, paramLocal);

%% options for fmincon
nParamLocal = length(optionsOptimLocal.initParams);
  
nMaxEvalLocalOpt = max([...
    optionsOptimLocal.nMaxEvalLocalOpt, ...
    ceil(optionsMFA.optionsOptimMH.nMaxEval/10)]);
optionsFmincon = optimset('Display', 'Iter',...
    'MaxIter', 10^4, 'MaxFunEvals', nMaxEvalLocalOpt,...
    'TolX', 10^-9);
nonlconFmincon = nonlcon;

%% Objective function
objfunFmincon = objfun;

%% Run local optimization
disp([' Maximum number of evaluation: ' num2str(nMaxEvalLocalOpt)])

[solLocalRaw.params, solLocalRaw.score, ...
    solLocalRaw.exitFlag, solLocalRaw.output, ~, solLocalRaw.grad, solLocalRaw.hessian]  ...
    = fmincon(objfunFmincon, initParamLocal, ...
    A, b, Aeq, beq, ...
    optionsOptimLocal.lb, optionsOptimLocal.ub, nonlconFmincon, optionsFmincon);


disp(' ')
disp('Calculating jacobian..')
optionsCalcJacobian.evalFxnType = 'obj';
try
    optionsCalcJacobian.diffType = {'central', 'forward'};
    tmpSol.jacobian = calcJacobian(objfunRes, solLocalRaw.params, ...
        optionsOptimLocal.lb, optionsOptimLocal.ub, optionsCalcJacobian);
%     save(optionsMFA.saveFileName,'optionsCalcJacobian','-append')
catch err
    tmpSol.jacobian = [];
end
if ~isempty(tmpSol.jacobian)
    disp('Jacobian was calculated.')
end

%% Final solution of local optimization
[~,~,solLocal] = objfun(solLocalRaw.params);
solLocal.jacobian = tmpSol.jacobian;
optionsOptimLocal.optionsLsqnonlinFmincon = optionsFmincon;
solLocal.output = solLocalRaw.output;

end


%% Function to calculate residuals
function residual=calcRes(objfun, params)
[~,~,~,residual] = objfun(params);
end

%% Make QP model for linear constrains for fmincon
function [qpModel] = makeQpModelFmincon(model, expData, optionsMFA, initParamLocal)

optimType = 'local';

isInputIndVars = true;
isInputFullIndVarsConstr= true;
fxnPrepQpModel = str2func(optionsMFA.fxnName.prepQpModel);
[qpModel] = fxnPrepQpModel(model, expData, optionsMFA,  ...
    optimType, initParamLocal,isInputIndVars, isInputFullIndVarsConstr);

end

%% Non-linear constrains for fmincon
function [C, Ceq] = makeNonLinConstrFmincon(model, expData, optionsMFA, paramLocal)

optimType = 'localNonLinConstr';

isInputIndVars = true;
isInputFullIndVarsConstr= true;
fxnPrepQpModel = str2func(optionsMFA.fxnName.prepQpModel);
[qpModel] = fxnPrepQpModel(model, expData, optionsMFA,  ...
    optimType, paramLocal,isInputIndVars, isInputFullIndVarsConstr);

C = qpModel.A * paramLocal - qpModel.b;
Ceq = [];

end



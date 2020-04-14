function [sol, model, expData, optionsMFA] =  MFA_NSS(model, expData, optionsMFA)
%MFA_NSS        run metabolic flux analysis under non-steady-state
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

import packFxnMFANSS.*

%% Load default options
defaultOptionsMFA = loadDefaultOptionsMFA();
defaultOptionsOptimMH = defaultOptionsMFA.optionsOptimMH;
filedNamesInputOptionsMFA = fieldnames(optionsMFA);
for f = 1 : length(filedNamesInputOptionsMFA)
    defaultOptionsMFA.(filedNamesInputOptionsMFA{f}) = ...
        optionsMFA.(filedNamesInputOptionsMFA{f});
end
filedNamesInputOptionsOptimMH = fieldnames(optionsMFA.optionsOptimMH);
for f = 1 : length(filedNamesInputOptionsOptimMH)
    defaultOptionsOptimMH.(filedNamesInputOptionsOptimMH{f}) = ...
        optionsMFA.optionsOptimMH.(filedNamesInputOptionsOptimMH{f});
end
optionsMFA = defaultOptionsMFA;
optionsMFA.optionsOptimMH = defaultOptionsOptimMH;

if optionsMFA.isRxnDepNSwitchTimes
    optionsMFA.nSwitchTimes = max(optionsMFA.nSwitchTimesNetRxns);
    optionsMFA.patternNSwitchTimes = unique(optionsMFA.nSwitchTimesNetRxns);
else
    optionsMFA.patternNSwitchTimes = optionsMFA.nSwitchTimes;
end

%% Add information of model
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
[idIndVars, idRedandantConstr, Cx, C0, idInvalidIndFluxes] = ...
    selectIndVars(A, isInputFullIndVarsConstr, idInputIndVars);
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

optionsMFA.transformMat = prepTransformMat(model, expData, optionsMFA);

expData.vecExpConcStruct = makeVecExpConcs(model, expData, optionsMFA);

%% Identify independent parameters
[optionsMFA.isIndParams.init, optionsMFA.isRedundantConstrAeq.init] = ...
    identifyIndMFAParams(model, expData, optionsMFA, 'init');
[optionsMFA.isIndParams.MH, optionsMFA.isRedundantConstrAeq.MH] = ...
    identifyIndMFAParams(model, expData, optionsMFA, 'metaheuristic');
[optionsMFA.isIndParams.local, optionsMFA.isRedundantConstrAeq.local] = ...
    identifyIndMFAParams(model, expData, optionsMFA, 'local');

optionsMFA.transformMat = prepTransformMat(model, expData, optionsMFA);
if optionsMFA.isUseQPInMH 
    [optionsMFA.isIndParams.QP, optionsMFA.isRedundantConstrAeq.QP] = ...
        identifyIndMFAParams(model, expData, optionsMFA, 'QP');
end

%% Options for ODE 
for i = 1 : length(optionsMFA.optionsODE.activeOptionNames)
    tmpOptionName = optionsMFA.optionsODE.activeOptionNames{i};
    if i == 1
        optionsMFA.optionsODE.final  = odeset(tmpOptionName, optionsMFA.optionsODE.(tmpOptionName));
    else
        optionsMFA.optionsODE.final= odeset(optionsMFA.optionsODE.final, tmpOptionName, optionsMFA.optionsODE.(tmpOptionName));
    end
end

disp('MFA-NSS start')
if optionsMFA.isRxnDepNSwitchTimes
    disp(['reaction-specific # of time intervals'])
else
    disp(['# of time intervals: ' num2str(optionsMFA.nSwitchTimes+1) ...
        ' (' num2str(optionsMFA.nSwitchTimes) ' switch times)'])
end
disp(['# of total parameters: ' num2str(optionsMFA.idParamLocal.nParam)])
disp(['# of parameters in metaheuristic optimization: ' num2str(optionsMFA.idParamMH.nParam)])
disp(' ')

%% Initial parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[initParams, initScores, nExpData] =  makeInitParams(...
    model, expData,optionsMFA);

expData.nExpData = nExpData;
optionsMFA.optionsOptimMH.initParams = initParams;
optionsMFA.optionsOptimMH.initScores = initScores;

disp(['Best Score by initial sampling before optimization:' num2str(min(initScores))])




%% Metaheuristic optimization %%%%%%%%%%%%%%%%%%%%%%%%%%%%
optimType =  'metaheuristic';
objfun =@(params) solveMDVDynamics(...
    params, model, expData, optionsMFA, optimType);
disp(' ')
disp('Metaheuristic optimization by CMA-ES')
solMH = performOptimMH(model, expData, optionsMFA, objfun, initParams);
disp(['Best Score by metaheuristic optimization:' num2str(solMH.score)])

%% Local optimization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
optimType =  'local';
objfun =@(params) solveMDVDynamics(...
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
optionsOptimMH.lb = param2vector(lbMH, optionsMFA, optimType, isLogTransform, isIndParam);
optionsOptimMH.ub = param2vector(ubMH, optionsMFA, optimType, isLogTransform, isIndParam);

%% Setting initial points for CMA-ES
% mean
meanInitParams = mean(initParams, 2);

% covariance matrix
tmpCov = cov(initParams');
diagCov = diag(tmpCov);
tmpSigma = sqrt(diagCov);

[V,D] = eig(tmpCov);
diagD = diag(D);
epsilon = 1e-12;
diagD(diagD<=epsilon) = epsilon;
D = diag(diagD);
tmpCov_old = tmpCov;
tmpCov = V * D * inv(V);

%% Setting options of CMA-ES
objfunStr = 'solveMDVDynamics';
if optionsOptimMH.isUseInitCovMat
    sigmaCMAES = tmpCov;
else
    sigmaCMAES = tmpSigma;
end
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
    cmaes_modified(objfunStr, meanInitParams, sigmaCMAES, optionsCMAES, varaginCMAES{:});

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

optionsOptimLocal.lb = param2vector(lbLocal, optionsMFA, optimType, isLogTransform, isIndParam);
optionsOptimLocal.ub = param2vector(ubLocal, optionsMFA, optimType, isLogTransform, isIndParam);

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

nMaxEvalLocalOpt = optionsOptimLocal.nMaxEvalLocalOpt;
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
[qpModel] = prepQpModel(model, expData, optionsMFA,  ...
    optimType, initParamLocal,isInputIndVars, isInputFullIndVarsConstr);

end

%% Non-linear constrains for fmincon
function [C, Ceq] = makeNonLinConstrFmincon(model, expData, optionsMFA, paramLocal)

optimType = 'localNonLinConstr';

isInputIndVars = true;
isInputFullIndVarsConstr= true;
[qpModel] = prepQpModel(model, expData, optionsMFA,  ...
    optimType, paramLocal,isInputIndVars, isInputFullIndVarsConstr);

C = qpModel.A * paramLocal - qpModel.b;
Ceq = [];

end



%% Run identification of effective flux regulators
function [solFinal, solInit, optionsIEFR] = run_IEFR(expDataKA, optionsIEFR)
%run_IEFR  run identification of effective flux regulators
% INPUT
% expDataKA     data structure of flux, substrates, products, phosphorylation, and allsoteric effectors
%        regulatorType      type of regulator
%                                 either of flux, substrates, products, phospho (phosphorylation), 
%                                 alloA (allosteric activation), or alloI (allosteric inhibition)
%        regulatorName     regulator name (e.g., glucose, Gys_S728)
%        ScoefRegulatorMet   stoichiometric coefficient (applied to only substrates and products)
%        level       amount of regulator
% optionsIEFR     options for identification of effective flux regulator
%
% OUTPUT
% solFinal      solution structure after minimization of  AIC
%       param       estimated parameters
%       paramInfo       information of the paramters
%       residuals       residuals between input and estimated fluxes
%       fluxEst         estimated fluxes
%       nExpData        number of input data point
%       nParam          number of parameters
%       AIC         Akaike's information criterion
%       effectiveFluxRegulators         identified effective flux regulators of phosphorylation and allosteric effectors
%       effectiveFluxRegulatorTypes         type of identified effective flux regulators
% solInit      solution structure before minimization of AIC
% optionsIEFR     options where default options, "normExpData" and "paramInfo" are added
%
% see example_IEFR.m as an example code.


import packFxnKARAMCAInit.*

%% Load default options
defaultOptionsKA = loadDefaultOptionsIEFR();
filedNamesInputOptionsKA= fieldnames(optionsIEFR);
for f = 1 : length(filedNamesInputOptionsKA)
    defaultOptionsKA.(filedNamesInputOptionsKA{f}) = optionsIEFR.(filedNamesInputOptionsKA{f});
end
optionsIEFR = defaultOptionsKA;

%% Normalize data by the L2 norm
normExpData= makeNormExpData(expDataKA, optionsIEFR);
optionsIEFR.normExpData = normExpData;

%% Preparation of parameter information
paramInfo.paramTypes = {'k1', 'k2', 'Ks', 'Kp', 'Kr'}';
paramInfo.paramRegulatorTypes = {'constant',  'phospho', 'alloA', 'alloI', 'substrate', 'product'}';
paramInfo.fxnGetOriParam = cell(length(paramInfo.paramRegulatorTypes),1);
for p = 1 : length(paramInfo.paramRegulatorTypes)
    switch paramInfo.paramRegulatorTypes{p}
        case {'phospho'}
            paramInfo.fxnGetOriParam{p} = '@(x) x';
        case {'alloI'}
            paramInfo.fxnGetOriParam{p} = '@(x) -10^x';
        otherwise
            paramInfo.fxnGetOriParam{p} = '@(x) 10^x';
    end
end
paramInfo.paramRegulatorNames = {'kcat+', 'kcat-', normExpData(2:end).regulatorName}';
optionsIEFR.paramInfo = paramInfo;

%% Prepalation of stepwise model selection
loc = ~ismember({normExpData.regulatorType}, {'flux'});
Xall = vertcat(normExpData(loc).level)';
Xall = [ones(size(Xall,1),2), Xall];
regulatorTypes = [{'constant','constant'},{normExpData(2:end).regulatorType}];
nParam = size(Xall, 2);
y = normExpData(1).level';
objfunIter = @(isUseVar) calcScoreIter(...
    y, Xall, isUseVar, normExpData, optionsIEFR);
objfunStepwiseFit = @(isUserVar) calcScoreStepwiseFit(...
    objfunIter, isUserVar, optionsIEFR);

%% Set which variable can be selected or not
% If both lb(i) and ub(i) equal to zero, variable i is never selected. 
% If both lb(i) and ub(i) equal to one, variable i is always selected.
% If both lb(i) equals to zero and ub(i) equals to one, variable i can be either selected or not.

options.lb = zeros(nParam, 1);
options.ub = ones(nParam, 1);

isAlwaysAddVar = false(nParam,1);

%Substrate and products are always selected.
isSubsProd = ismember(regulatorTypes, {'substrate', 'product'});
isAlwaysAddVar(isSubsProd) = 1;

% k1 is always selected.
isAlwaysAddVar(1) = true; 
% k2 is always selcted only if product exists
if any(strcmp(regulatorTypes, 'product'))
    isAlwaysAddVar(2) = true;
else
    options.ub(2) = 0;
end
options.lb(isAlwaysAddVar) = 1;


%% Run stepwise model selection
[solStepwiseFinal, solStepwiseInit] = myStepwiseSelection(Xall, objfunStepwiseFit, options);

%% Initial model
solInit = solStepwiseInit.sol;
solInit = rmfield(solInit, {'jacobian', 'optimLocal', 'isUseVar', 'optionsOptimMH', 'optionsOptimLocal', 'y', 'X', 'score'});
solInit.paramInfo.paramRegulatorTypes = ...
    optionsIEFR.paramInfo.paramRegulatorTypes(solInit.paramInfo.paramRegulatorTypeIds);
solInit.paramInfo.paramRegulatorNames = ...
    optionsIEFR.paramInfo.paramRegulatorNames(solInit.paramInfo.paramRegulatorNameIds);

idEFRNormExpData = solInit.paramInfo.locNormExpData(solInit.paramInfo.isRegParam);
if isempty(idEFRNormExpData)
    solInit.effectiveFluxRegulators = {};
    solInit.effectiveFluxRegulatorTypes = {};
else
    solInit.effectiveFluxRegulators = {optionsIEFR.normExpData(idEFRNormExpData).regulatorName}';
    solInit.effectiveFluxRegulatorTypes = {optionsIEFR.normExpData(idEFRNormExpData).regulatorType}';
end


%% Final model
solFinal = solStepwiseFinal.sol;
solFinal = rmfield(solFinal, {'jacobian', 'optimLocal', 'isUseVar', 'optionsOptimMH', 'optionsOptimLocal', 'y', 'X', 'score'});
solFinal.paramInfo.paramRegulatorTypes = ...
    optionsIEFR.paramInfo.paramRegulatorTypes(solFinal.paramInfo.paramRegulatorTypeIds);
solFinal.paramInfo.paramRegulatorNames = ...
    optionsIEFR.paramInfo.paramRegulatorNames(solFinal.paramInfo.paramRegulatorNameIds);

idEFRNormExpData = solFinal.paramInfo.locNormExpData(solFinal.paramInfo.isRegParam);
if isempty(idEFRNormExpData)
    solFinal.effectiveFluxRegulators = {};
    solFinal.effectiveFluxRegulatorTypes = {};
else
    solFinal.effectiveFluxRegulators = {optionsIEFR.normExpData(idEFRNormExpData).regulatorName}';
    solFinal.effectiveFluxRegulatorTypes = {optionsIEFR.normExpData(idEFRNormExpData).regulatorType}';
end

end

%% Run independent parameter estimation
function [score, sol]= calcScoreStepwiseFit(objfun, isUseVar, optionsIEFR)

nIterOpt = optionsIEFR.nIterOpt;
scoresAllIter = zeros(1, nIterOpt);
solIter = cell(1, nIterOpt);
for i = 1 : nIterOpt
    [scoresAllIter(i), solIter{i}] = objfun(isUseVar);
end

[~, idMin] = min(scoresAllIter);
score = scoresAllIter(idMin);
sol = solIter{idMin};

end

%% Calculation of score at each parameter estimation
function [score, sol] = calcScoreIter(y, X, isUseVar, normExpData, optionsIEFR)

regulatorNames = [{'kcat+','kcat-'},{normExpData(2:end).regulatorName}];
regulatorTypes = [{'constant','constant'},{normExpData(2:end).regulatorType}];

nParam = nnz(isUseVar);

%% Information of parameters

paramInfoAll = optionsIEFR.paramInfo;

paramTypeIds = zeros(nParam, 1);
paramRegulatorTypeIds = zeros(nParam, 1);
paramRegulatorNameIds = zeros(nParam, 1);
locNormExpData = zeros(nParam,1);
vv = 0;
for v = 1 : length(isUseVar)
    if ~isUseVar(v)
        continue
    end
    vv = vv + 1;
    switch regulatorTypes{v}
        case {'constant'} 
            if v==1% k1
                tmpParamTypeId = 1;
            else% k2
                tmpParamTypeId = 2;
            end
        case {'substrate'}
            tmpParamTypeId = 3;
        case {'product'}
            tmpParamTypeId = 4;
        otherwise % allosteric or phosphorylation
            tmpParamTypeId = 5;
    end
    tmpRegulatorType = regulatorTypes(v);
    tmpRegulatorName = regulatorNames(v);
    
    paramTypeIds(vv) = tmpParamTypeId;
    [TF, paramRegulatorTypeIds(vv)] = ismember(tmpRegulatorType, paramInfoAll.paramRegulatorTypes);
    if ~all(TF)
        error('paramRegulatorTypes is not correct.')
    end
    [TF, paramRegulatorNameIds(vv)] = ismember(tmpRegulatorName, paramInfoAll.paramRegulatorNames);
    if ~all(TF)
        error('paramRegulatorNames fis not correct.')
    end
    
    switch v
        case 1
            locNormExpData(vv) = v;
        case 2
            locNormExpData(vv) = 0;
        otherwise
            locNormExpData(vv) = v-1;
    end
end

paramInfo.paramTypeIds = paramTypeIds;
paramInfo.paramRegulatorTypeIds = paramRegulatorTypeIds;
paramInfo.paramRegulatorNameIds = paramRegulatorNameIds;
paramInfo.locNormExpData = locNormExpData;
[~, regParamRegulatorTypeIds] = ismember({'phospho','alloA', 'alloI'}, paramInfoAll.paramRegulatorTypes);
paramInfo.isRegParam = ismember(paramInfo.paramRegulatorTypeIds, regParamRegulatorTypeIds);


%% Set bound of parameters
optionsOpt = optionsIEFR.optionsOpt;
optionsOptimMH.popSize = optionsOpt.popSize;
optionsOptimMH.nMaxEval = optionsOpt.nMaxEval;
maxAbsVal = optionsOpt.maxAbsVal;
minAbsVal = optionsOpt.minAbsVal;

optionsOptimMH.ub = zeros(nParam,1);
optionsOptimMH.lb = zeros(nParam,1);
for p = 1 : nParam
    isKr = false;
    switch paramInfoAll.paramTypes{paramInfo.paramTypeIds(p)}
        case {'k1','k2'}
            tmpUb = log10(maxAbsVal)+10;
            tmpLb = log10(minAbsVal)-2;
        case {'Ks', 'Kp'}
            tmpUb = log10(maxAbsVal);
            tmpLb = log10(minAbsVal);
        case {'Kr'}
            isKr = true;
        otherwise
            error
    end
    
    if isKr
        switch paramInfoAll.paramRegulatorTypes{paramInfo.paramRegulatorTypeIds(p)}
            case {'phospho'}
                tmpUb = +maxAbsVal;
                tmpLb = -maxAbsVal;
            case {'alloA'}
                tmpUb = log10(maxAbsVal);
                tmpLb = log10(minAbsVal);
            case {'alloI'}
                tmpUb = log10(maxAbsVal);
                tmpLb = log10(minAbsVal);
            otherwise
                error
        end
    end
    
    optionsOptimMH.ub(p) = tmpUb;
    optionsOptimMH.lb(p) = tmpLb;
end

%% Score function
for p = 1 : length(optionsIEFR.paramInfo.fxnGetOriParam)
    optionsIEFR.paramInfo.fxnGetOriParam{p} = str2func(optionsIEFR.paramInfo.fxnGetOriParam{p});
end

fxnCalcRes = @(param) calcRes(y, X(:,isUseVar), param, normExpData, paramInfo, optionsIEFR);
objfunRSS =@(param) calcRSS(param, fxnCalcRes, optionsOpt.infeasibleScore);
objfunStr = 'calcRSS';

%% Initial parameters
idIndv = 0;
nFailure = 0;
isInfeasible=false;
optionsOptimMH.initScores = zeros(1,optionsOpt.popSize);
optionsOptimMH.initParams = zeros(nParam, optionsOpt.popSize);
while true
    idIndv = idIndv + 1;
    locParamPhospho = ismember(...
        paramInfoAll.paramRegulatorTypes(paramInfo.paramRegulatorTypeIds), {'phospho'});
    
    ubInitParams = optionsOptimMH.ub;
    lbInitParams = optionsOptimMH.lb;
    ubInitParams(locParamPhospho) =  log10(maxAbsVal);
    lbInitParams(locParamPhospho) = log10(minAbsVal);
    
    randVal = rand(nParam,1);
    paramRange = ubInitParams - lbInitParams;
    tmpParam = randVal.*paramRange + lbInitParams;
    if any(locParamPhospho)
        tmpParam(locParamPhospho) = 10.^ tmpParam(locParamPhospho);
        tmpSign = ones(nnz(locParamPhospho),1);
        tmpSign(rand(nnz(locParamPhospho),1)<0.5) = -1;
        tmpParam(locParamPhospho) = tmpParam(locParamPhospho) .* tmpSign;
    end
    
    tmpScore = objfunRSS(tmpParam);
    if isreal(tmpScore) && tmpScore<10^6
        optionsOptimMH.initParams(:,idIndv) = tmpParam;
        optionsOptimMH.initScores(idIndv) = tmpScore;
    else
        nFailure = nFailure + 1;
        idIndv=idIndv-1;
    end
    
    if idIndv>=optionsOpt.popSize
        break
    end
    if nFailure >= max([optionsOpt.nMaxEval, 10^4])
        isInfeasible = true;
        break
    end
    
end


%% Feasible or not
if isInfeasible
    warning(['Feasible initical points were not identified from ' num2str(nFailure) ' random parameter sets'])
    sol.param = [];
    sol.score = optionsOpt.infeasibleScore;
    sol.jacobian = [];
    sol.optimLocal = false;
    sol.paramInfo = [];
    sol.residuals = [];
    sol.fluxEst= [];
    sol.factor = [];
    sol.nExpData = nnz(y);
    sol.nParam = nParam;
    sol.AIC = sol.nExpData*log(sol.score/sol.nExpData) +  2 *nParam;
    sol.historyMH = [];
    sol.y = y;
    sol.X = X;
    sol.isUseVar = isUseVar;
    sol.optionsOptimMH = optionsOptimMH;
    sol.optionsOptimLocal = [];
    score = sol.AIC;
    return
end



%% Run CMA-ES
% saveFun = [];
sigmaCMAES = std(optionsOptimMH.initParams,[],2);
varaginCMAES{1} = fxnCalcRes;
varaginCMAES{2} = optionsOpt.infeasibleScore;
optionsCMAES = cmaes('defaults');
optionsCMAES.MaxFunEvals = optionsOptimMH.nMaxEval;
optionsCMAES.LBounds = optionsOptimMH.lb;
optionsCMAES.UBounds = optionsOptimMH.ub;
optionsCMAES.PopSize = optionsOptimMH.popSize;
optionsCMAES.DispModulo =inf;
optionsCMAES.DispFinal ='off';
optionsCMAES.CMA.active = 1;
optionsCMAES.SaveVariables = 'off';
if nParam >=2
    warning off
    [solCMAES.bestParam, solCMAES.bestScore, ...
        solCMAES.couteval, solCMAES.stopflag,solCMAES.out,solCMAES.bestever] = ...
        cmaes(objfunStr, optionsOptimMH.initParams, sigmaCMAES, optionsCMAES, varaginCMAES{:});
    warning on
    sol.param = solCMAES.bestever.x;
    sol.score = solCMAES.bestever.f;
    sol.jacobian = [];
    solMH = sol;
    clear sol
else
    scoreAll = zeros(1,optionsOptimMH.popSize);
    for p = 1 : optionsOptimMH.popSize
        scoreAll(p) = objfunRSS(optionsOptimMH.initParams(p));
    end
    [~, idMin] = min(scoreAll);
    sol.param = optionsOptimMH.initParams(idMin);
    sol.score = scoreAll(idMin);
    sol.jacobian = [];
    solMH = sol;
    historyMH = [];
    clear sol
end

%% Run fmincon
optionsOptimLocal = optimset('Display', 'none');
optionsOptimLocal.initParams = solMH.param;
optionsOptimLocal.initScores = solMH.score;
optionsOptimLocal.lb = optionsOptimMH.lb;
optionsOptimLocal.ub = optionsOptimMH.ub;

locParamPhospho = ismember(...
    paramInfoAll.paramRegulatorTypes(paramInfo.paramRegulatorTypeIds), {'phospho'});
% phosphorylation
optionsOptimLocal.lb(locParamPhospho) = optionsOptimLocal.lb(locParamPhospho)*10;
optionsOptimLocal.ub(locParamPhospho) = optionsOptimLocal.ub(locParamPhospho)*10;
% others
optionsOptimLocal.lb(~locParamPhospho) = optionsOptimLocal.lb(~locParamPhospho)-1;
optionsOptimLocal.ub(~locParamPhospho) = optionsOptimLocal.ub(~locParamPhospho)+1;

[solLocal.param, solLocal.score, ~, ...
    ~, ~, ~, solLocal.jacobian]  ...
    = fmincon(objfunRSS, optionsOptimLocal.initParams, ...
    [],[],[],[],optionsOptimLocal.lb, optionsOptimLocal.ub, [],optionsOptimLocal);        

%% Update solution
if solLocal.score > solMH.score
    sol = solMH;
    sol.optimLocal = false;
else
    sol = solLocal;
    sol.optimLocal = true;
end
sol.paramInfo = paramInfo;

[~, sol.residuals, sol.fluxEst, sol.factor] = objfunRSS(sol.param);
if sol.factor.fR<0
    error('fR must be positive')
end

sol.nExpData = nnz(y);
sol.nParam = nParam;

sol.AIC = sol.nExpData*log(sol.score/sol.nExpData) +  2 *nParam;
sol.y = y;
sol.X = X;
sol.isUseVar = isUseVar;
sol.optionsOptimMH = optionsOptimMH;
sol.optionsOptimLocal = optionsOptimLocal;

score = sol.AIC;

end


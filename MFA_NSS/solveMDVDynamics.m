%% Simulate MDV
function [score, scoreMet, sol, residual] = ...
    solveMDVDynamics171206(paramVector, model, expData, optionsMFA, optimType, scoreInfeasible)
% function [score, scoreMet, MDVsSim, timeSim, concs, nExpData, isCompMat] = ...
%     solveMDVDynamics160421(model, expData, paramVector, optionsMFA)
% IDV*dC/dt + C*dIDV/dt = ƒ°S*v*IMM*IDV_i' - ƒ°(-S*v*IDV_i)
% MDV = ICM * IDV
% 
warning off
field2var(optionsMFA.varSet)
idParamLocal = optionsMFA.idParamLocal;
nParamLocal = idParamLocal.nParam;
if nargin<=5
    scoreInfeasible = 10^12;
end

%% Prepare return values
score = scoreInfeasible;
scoreMet = zeros(nEvalMDVMets,1) + score/nEvalMDVMets;
sol =[];
residual = [];
if isfield(expData, 'nExpData')
    residual = zeros(expData.nExpData,1) + sqrt(score/expData.nExpData);
end

%% Check difference between knot times
switch optimType
    case {'local'}
        idParam = optionsMFA.idParamLocal;
        isInParams = false(1,idParam.nParam);
        isInParams(idParam.knots) = true;
        isInParams = isInParams(optionsMFA.isIndParams.local);
        knots = paramVector(isInParams);
        
    case {'metaheuristic', 'init'}
        idParam = optionsMFA.idParamMH;
        isInParams = false(1,idParam.nParam);
        isInParams(idParam.knots) = true;
        isInParams = isInParams(optionsMFA.isIndParams.MH);
        knots = 10.^paramVector(isInParams);        
end

fullKnots = [0, knots', expData.time(end)];
if any(diff(fullKnots) < optionsMFA.minKnotTimeDiff)
    return
end

%% Prepare paramLocal
switch optimType
    case {'metaheuristic', 'init'}
            tmpOptimType = 'QP';
    otherwise
        tmpOptimType = optimType;
end
fxnParamInd2LocalFull = str2func(optionsMFA.fxnName.paramInd2LocalFull);
[paramLocal, convertMat] = fxnParamInd2LocalFull(model, expData, optionsMFA, paramVector, tmpOptimType);
if isempty(paramLocal)
    return
end
optionsMFA.convertMat = convertMat;

%% Flux
fxnVector2param = str2func(optionsMFA.fxnName.vector2param);
[knotFluxes, initConcs, knotConcRates, knots, coefCorrCompt] = ...
    fxnVector2param(paramLocal, model, expData, optionsMFA);
if isempty(knotFluxes)
    return
end

if any(any(knotFluxes < optionsMFA.lbLocal.knotFluxesAll*(1-10^-3)))|| ...
        any(any(knotFluxes > optionsMFA.ubLocal.knotFluxesAll*(1+10^-3)))
    return
end


%% Metabolite concentrations
concs = convertMat.param2Conc * paramLocal; 

nTmpConstr = size(convertMat.param2Conc, 1);
tmpLb = diag(optionsMFA.lbLocal.initConcs)*ones(nNonPoolMets, optionsMFA.nStepTimeConstrQP);

tmpLb = tmpLb/optionsMFA.paramBoundExpandRatio;
if any(concs < reshape(tmpLb', nTmpConstr, 1)*(1-10^-3))
    return
end

concsAns = reshape(concs, optionsMFA.nStepTimeConstrQP, nNonPoolMets)';
initConcsTest = concsAns(:,1);
initConcsAns = convertMat.param2InitConc * paramLocal;

%% polynomial coefficient of metabolite concetnrations
pCoefConcs = convertMat.param2pCoefConc*paramLocal;
pCoefConcs = reshape(pCoefConcs, nNonPoolMets, 3*(nKnots+1));


%% polynomial coefficient of flux

pCoefFluxes = convertMat.param2pCoefFlux * paramLocal;
pCoefFluxes = reshape(pCoefFluxes,nRxns, 2*(nKnots+1));


%% Prepare EMU simulation
for k = 1 : nKnots + 1
    if nKnots == 0
        timeSpan = fullKnots;
    else
        timeSpan = [fullKnots(k), fullKnots(k+1)];
    end
    isExpTimeInTimeSpan = expData.time>=timeSpan(1) & expData.time<=timeSpan(end);
    timeSpan = sort(unique([timeSpan, expData.time(isExpTimeInTimeSpan)]));
    
    if timeSpan(1) ==0
        if length(timeSpan) == 2
            querryTime = linspace(timeSpan(1), timeSpan(2),200);
        else
            querryTime = [linspace(timeSpan(1), timeSpan(2),200), linspace(timeSpan(2)*1.001, timeSpan(end),100)];
        end
    else
        querryTime = [linspace(timeSpan(1), timeSpan(end),100)];
    end
    if length(timeSpan)>=3
        querryTime = sort(unique([querryTime, timeSpan(2:end-1)]));
    end
    
    tmpParamStruct.knotFluxes = knotFluxes(:,k+(0:1));
    tmpParamStruct.knots = knots;
    tmpParamStruct.pCoefConcs = pCoefConcs(:, (k-1)*3+(1:3));
    tmpParamStruct.pCoefFluxes = pCoefFluxes(:, (k-1)*2+(1:2));
    tmpParamStruct.querryTime= querryTime;
    
    concsQuerryTime = tmpParamStruct.pCoefConcs*[querryTime.^0; querryTime.^1; querryTime.^2];
    if any(any(concsQuerryTime < repmat(optionsMFA.lbLocal.initConcs/10, 1, length(querryTime))))
        return
    end
    paramStruct(k) = tmpParamStruct;
end



%% Simulate EMU
nTimeSim = zeros(1,nKnots+1);

fxnSoveEMUDynamics = str2func(optionsMFA.fxnName.solveEMUDynamics);
for k = 1 : nKnots + 1
    try
        if k == 1
            output(k) = fxnSoveEMUDynamics(paramStruct(k), model, expData, optionsMFA, []);
        else
            output(k) = fxnSoveEMUDynamics(paramStruct(k), model, expData, optionsMFA, output(k-1));
        end
        if isempty(output(k))
            return
        end
        nTimeSim(k) = length(output(k).timeSim);
        if isempty(output(k).concs)
            return
        end
        if any(output(k).concs(:,end)<0) || output(k).timeSim(end) < paramStruct(k).querryTime(end)
            return
        end
    catch err
        return
    end
end


%% merge results

MDVsSim = cell(nNonPoolMets, 1);
for m = 1 : nMets
        tmpMDVsSim = cell(1,nKnots+1);
        for k = 1:nKnots+1
            tmpMDVsSim(k) = output(k).MDVsSim(m);
        end
        MDVsSim{m} = [tmpMDVsSim{:}];
        MDVsSim{m} = single(MDVsSim{m});
end
clear tmpMDVsSim

% Correct natural abundance
if optionsMFA.isCorrectNatAbundanceSim
    for m = 1 : nNonPoolMets
        tmpMDVsSim = MDVsSim{m};
        nMDV = size(tmpMDVsSim,1);
        corrMDVsSim = expData.corrMat{m}(1:nMDV,1:nMDV) * tmpMDVsSim;
        sumMDVsSim = sum(corrMDVsSim ,1);
        corrMDVsSim = corrMDVsSim ./ repmat(sumMDVsSim, nMDV, 1);
        MDVsSim{m} = corrMDVsSim;
        MDVsSim{m} = single(MDVsSim{m});
    end
end

if size(output(1).timeSim,1)>size(output(1).timeSim,2)
    timeSim = vertcat(output.timeSim)';
else
    timeSim = [output.timeSim];
end
concs = [output.concs];

concMDVs = cell(size(MDVsSim));
for m = 1 : nNonPoolMets
        concMDVs{m} = MDVsSim{m} .* repmat(concs(m,:),size(MDVsSim{m},1),1);
        concMDVs{m} = single(concMDVs{m});
end

sol.knotFluxes = knotFluxes;
sol.netFluxes = model.matRaw2Net*knotFluxes; 
sol.knots = knots;
sol.MDVsSim = MDVsSim;
sol.concs = single(concs);
sol.concMDVs = concMDVs;
sol.initConcs = initConcs;
sol.knotConcRates = knotConcRates;
sol.pCoefConcs = pCoefConcs;
sol.timeSim = timeSim;

if ~isempty(idComptCorrParam)
    sol.coefCorrCompt = coefCorrCompt;
else
        mediaInfo = optionsMFA.mediaInfo;
        sol.coefCorrCompt = 1/mediaInfo.mediaVol*mediaInfo.cellAmount;
end

%% RSS
[score, scoreMet, residual, resVec2MDV]= calcScore(model, expData, optionsMFA, sol);

sol.score = score;
sol.paramVector = paramLocal;
sol.residual = residual;
sol.resVec2MDV = resVec2MDV;
sol.nExpData = length(residual);

% sort
tmpFieldNames = fieldnames(sol);
firstFieldNames = {'score'; 'paramVector'};
tmpFieldNames = [firstFieldNames; tmpFieldNames(~ismember(tmpFieldNames, firstFieldNames))];
sol = orderfields(sol, tmpFieldNames);
warning on

end


%% Simulate MDV
function [score, scoreMet, sol, residual] = ...
    solveMDVDynamics(paramVector, model, expData, optionsMFA, optimType)

warning off
field2var(optionsMFA.varSet)
idParamLocal = optionsMFA.idParamLocal;
nParamLocal = idParamLocal.nParam;
if nargin<=5
    scoreInfeasible = nan;
end

%% Prepare return values
score = scoreInfeasible;
scoreMet = zeros(nEvalMDVMets,1) + score/nEvalMDVMets;
sol =[];
residual = [];
if isfield(expData, 'nExpData')
    residual = zeros(expData.nExpData,1) + sqrt(score/expData.nExpData);
end

%% Check difference between switch times
switch optimType
    case {'local'}
        idParam = optionsMFA.idParamLocal;
        isInParams = false(1,idParam.nParam);
        isInParams(idParam.switchTimes) = true;
        isInParams = isInParams(optionsMFA.isIndParams.local);
        switchTimes = paramVector(isInParams);
        
    case {'metaheuristic', 'init'}
        idParam = optionsMFA.idParamMH;
        isInParams = false(1,idParam.nParam);
        isInParams(idParam.switchTimes) = true;
        isInParams = isInParams(optionsMFA.isIndParams.MH);
        switchTimes = 10.^paramVector(isInParams);        
end

fullSwitchTimes = [0, switchTimes', expData.time(end)];
if any(diff(fullSwitchTimes) < optionsMFA.minKnotTimeDiff)
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
[switchTimeFluxes, initConcs, switchTimeConcRates, switchTimes] = ...
    fxnVector2param(paramLocal, model, expData, optionsMFA);
if isempty(switchTimeFluxes)
    return
end

if any(any(switchTimeFluxes < optionsMFA.lbLocal.switchTimeFluxesAll*(1-10^-3)))|| ...
        any(any(switchTimeFluxes > optionsMFA.ubLocal.switchTimeFluxesAll*(1+10^-3)))
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
pCoefConcs = reshape(pCoefConcs, nNonPoolMets, 3*(nSwitchTimes+1));


%% polynomial coefficient of flux

pCoefFluxes = convertMat.param2pCoefFlux * paramLocal;
pCoefFluxes = reshape(pCoefFluxes,nRxns, 2*(nSwitchTimes+1));


%% Prepare EMU simulation
for k = 1 : nSwitchTimes + 1
    if nSwitchTimes == 0
        timeSpan = fullSwitchTimes;
    else
        timeSpan = [fullSwitchTimes(k), fullSwitchTimes(k+1)];
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
    
    tmpParamStruct.switchTimeFluxes = switchTimeFluxes(:,k+(0:1));
    tmpParamStruct.switchTimes = switchTimes;
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
nTimeSim = zeros(1,nSwitchTimes+1);

fxnSoveEMUDynamics = str2func(optionsMFA.fxnName.solveEMUDynamics);
for k = 1 : nSwitchTimes + 1
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
        tmpMDVsSim = cell(1,nSwitchTimes+1);
        for k = 1:nSwitchTimes+1
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

% concMDVs = cell(size(MDVsSim));
% for m = 1 : nNonPoolMets
%         concMDVs{m} = MDVsSim{m} .* repmat(concs(m,:),size(MDVsSim{m},1),1);
%         concMDVs{m} = single(concMDVs{m});
% end

sol.switchTimeFluxes = switchTimeFluxes;
sol.netFluxes = model.matRaw2Net*switchTimeFluxes; 
sol.switchTimes = switchTimes;
sol.MDVsSim = MDVsSim;
sol.concs = single(concs);
% sol.concMDVs = concMDVs;
sol.initConcs = initConcs;
sol.switchTimeConcRates = switchTimeConcRates;
sol.pCoefConcs = pCoefConcs;
sol.timeSim = timeSim;

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


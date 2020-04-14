%% Simulate mass isotopomer fractions (MDV)
function [score, scoreMet, sol, residual] = ...
    solveMDVDynamics(paramVector, model, expData, optionsMFA, optimType)

warning off
nEvalMDVMets = optionsMFA.varSet.nEvalMDVMets;
nSwitchTimes = optionsMFA.varSet.nSwitchTimes;
nNonPoolMets = optionsMFA.varSet.nNonPoolMets;
nMets = optionsMFA.varSet.nMets;
nRxns = optionsMFA.varSet.nRxns;
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
[paramLocal, transformMat] = paramInd2LocalFull(model, expData, optionsMFA, paramVector, tmpOptimType);
if isempty(paramLocal)
    return
end
optionsMFA.transformMat = transformMat;

%% Flux
paramStruct = vector2paramStruct(paramLocal, model, expData, optionsMFA);
field2var(paramStruct)
if isempty(switchTimeFluxes)
    return
end

if any(any(switchTimeFluxes < optionsMFA.lbLocal.switchTimeFluxesAll*(1-10^-6)))|| ...
        any(any(switchTimeFluxes > optionsMFA.ubLocal.switchTimeFluxesAll*(1+10^-6)))
    return
end


%% Metabolite concentrations
concsQP = transformMat.param2ConcQP * paramLocal;
lbConcs = min(optionsMFA.lbLocal.switchTimeConcs, [],2);
tmpLb = diag(lbConcs)*ones(nNonPoolMets, optionsMFA.nStepTimeConstrQP);
tmpLb = tmpLb/optionsMFA.paramBoundExpandRatio;
tmpLb = tmpLb';
if any(concsQP < tmpLb(:)*(1-10^-6))
    return
end
ubConcs = min(optionsMFA.ubLocal.switchTimeConcs, [],2);
tmpUb = diag(ubConcs)*ones(nNonPoolMets, optionsMFA.nStepTimeConstrQP);
tmpUb = tmpUb*optionsMFA.paramBoundExpandRatio;
tmpUb = tmpUb';
if any(concsQP > tmpUb(:)*(1+10^-6))
    return
end

%% Concentration changes over time
if any(any(switchTimeConcRates < optionsMFA.lbLocal.switchTimeConcRates*(1-10^-6)))|| ...
        any(any(switchTimeConcRates > optionsMFA.ubLocal.switchTimeConcRates*(1+10^-6)))
    return
end

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
    
    tmpParamStructTI.switchTimeFluxes = switchTimeFluxes(:,k+(0:1));
    tmpParamStructTI.switchTimes = switchTimes;
    tmpParamStructTI.pCoefConcs = pCoefConcs(:, (k-1)*3+(1:3));
    tmpParamStructTI.pCoefFluxes = pCoefFluxes(:, (k-1)*2+(1:2));
    tmpParamStructTI.querryTime= querryTime;
    
    concsQuerryTime = tmpParamStructTI.pCoefConcs*[querryTime.^0; querryTime.^1; querryTime.^2];
    if any(any(concsQuerryTime < repmat(optionsMFA.lbLocal.initConcs/10, 1, length(querryTime))))
        return
    end
    paramStructTI(k) = tmpParamStructTI;
end



%% Simulate EMU
nTimeSim = zeros(1,nSwitchTimes+1);

for k = 1 : nSwitchTimes + 1
        if k == 1
            output(k) = solveEMUDynamics(paramStructTI(k), model, expData, optionsMFA, []);
        else
            output(k) = solveEMUDynamics(paramStructTI(k), model, expData, optionsMFA, output(k-1));
        end
        if isempty(output(k))
            return
        end
        nTimeSim(k) = length(output(k).timeSim);
        if isempty(output(k).concs)
            return
        end
        if any(output(k).concs(:,end)<0) || output(k).timeSim(end) < paramStructTI(k).querryTime(end)
            return
        end
end


%% Merge results

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

sol.switchTimeFluxes = switchTimeFluxes;
sol.netFluxes = model.matRaw2Net*switchTimeFluxes; 
sol.switchTimes = switchTimes;
sol.MDVsSim = MDVsSim;
sol.concs = single(concs);
sol.initConcs = paramStruct.initConcs;
sol.switchTimeConcs = paramStruct.switchTimeConcs;
sol.initConcRates = paramStruct.initConcs;
sol.switchTimeConcRates = switchTimeConcRates;
sol.pCoefConcs = paramStruct.pCoefConcs;
sol.pCoefFluxes= paramStruct.pCoefFluxes;
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


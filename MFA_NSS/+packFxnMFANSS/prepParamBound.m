%% Preparation of parameter bounds
function [lb, ub] = prepParamBound(model, expData, optionsMFA, optimType, expandRatio)

if nargin <= 4
    expandRatio = 1;
end

field2var(optionsMFA.varSet)

if ~optionsMFA.isUseQPInMH  
    if strcmp(optimType, 'metaheuristic')
        optimType = 'local';
    end
end

switch optimType
    case {'metaheuristic'}
        nOptimMetsConcs =  nMHMetsConcs;
        idOptimMetsConcs =idMHMetsConcs;
        nOptimMetsConcRates = nMHMetsConcRates;
        idOptimMetsConcRates = idMHMetsConcRates;
    case {'init', 'local'}
        nOptimMetsConcs = nNonPoolMets;
        idOptimMetsConcs =idNonPoolMets;
        nOptimMetsConcRates = nNonPoolMets;
        idOptimMetsConcRates =idNonPoolMets;
end


%% Bounds for initConcs, switchTimeConcRates, switchTimeFluxes, and switch times

metInfo = model.metInfo;
rxnInfo = model.rxnInfo;
lb.initConcs = metInfo.lbConcs(~metInfo.isPoolMets);
lb.switchTimeConcs = repmat(metInfo.lbConcs(~metInfo.isPoolMets),1, nSwitchTimes+2);
lb.initConcRates = metInfo.lbConcRates(~metInfo.isPoolMets);
lb.switchTimeConcRates = repmat(metInfo.lbConcRates(~metInfo.isPoolMets),1, nSwitchTimes+2);
lb.switchTimeFluxes = repmat(rxnInfo.lbFluxes(idIndFluxes),1, nSwitchTimes+2);
lb.switchTimes =  zeros(nSwitchTimes,1)+optionsMFA.lb.switchTimes+optionsMFA.minKnotTimeDiff;

ub.initConcs = metInfo.ubConcs(~metInfo.isPoolMets);
ub.switchTimeConcs = repmat(metInfo.ubConcs(~metInfo.isPoolMets),1, nSwitchTimes+2);
ub.initConcRates = metInfo.ubConcRates(~metInfo.isPoolMets);
ub.switchTimeConcRates = repmat(metInfo.ubConcRates(~metInfo.isPoolMets),1, nSwitchTimes+2);
ub.switchTimeFluxes = repmat(rxnInfo.ubFluxes(idIndFluxes),1, nSwitchTimes+2);
ub.switchTimes =  zeros(nSwitchTimes,1)+optionsMFA.ub.switchTimes-optionsMFA.minKnotTimeDiff;

switch optimType
    case {'metaheuristic'}
        lb.initConcs = lb.initConcs(idOptimMetsConcs,:);
        lb.switchTimeConcs = lb.switchTimeConcs(idOptimMetsConcs,:);
        lb.initConcRates = lb.initConcRates(idOptimMetsConcRates,:);
        lb.switchTimeConcRates = lb.switchTimeConcRates(idOptimMetsConcRates,:);
        ub.initConcs = ub.initConcs(idOptimMetsConcs,:);
        ub.switchTimeConcs = ub.switchTimeConcs(idOptimMetsConcs,:);
        ub.initConcRates = ub.initConcRates(idOptimMetsConcRates,:);
        ub.switchTimeConcRates = ub.switchTimeConcRates(idOptimMetsConcRates,:);
end


%% Expand parameter bound according to expandRatio
lb.initConcs = lb.initConcs/expandRatio;
lb.switchTimeConcs = lb.switchTimeConcs/expandRatio;
lb.initConcRates = lb.initConcRates*expandRatio;
lb.switchTimeConcRates= lb.switchTimeConcRates*expandRatio;
lb.switchTimeFluxes = lb.switchTimeFluxes/expandRatio;
lb.switchTimes =  lb.switchTimes -optionsMFA.minKnotTimeDiff+optionsMFA.minKnotTimeDiff/expandRatio;

ub.initConcs= ub.initConcs*expandRatio;
ub.switchTimeConcs= ub.switchTimeConcs*expandRatio;
ub.initConcRates= ub.initConcRates*expandRatio;
ub.switchTimeConcRates= ub.switchTimeConcRates*expandRatio;
ub.switchTimeFluxes = ub.switchTimeFluxes*expandRatio;
ub.switchTimes =  ub.switchTimes +optionsMFA.minKnotTimeDiff-optionsMFA.minKnotTimeDiff/expandRatio;

%% Bounds of fluxes
switch optimType
    case {'init','local'}
        lb.switchTimeFluxesAll =  repmat(rxnInfo.lbFluxes,1, nSwitchTimes+2);
        ub.switchTimeFluxesAll =  repmat(rxnInfo.ubFluxes,1, nSwitchTimes+2);
        lb.switchTimeFluxesAll = lb.switchTimeFluxesAll/expandRatio;
        ub.switchTimeFluxesAll = ub.switchTimeFluxesAll*expandRatio;
end


end
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
        nOptimMetsInitConcs =  nOuterOptimMetsInitConcs;
        idOptimMetsInitConcs =idOuterOptimMetsInitConcs;
        nOptimMetsConcRates = nOuterOptimMetsConcRates;
        idOptimMetsConcRates = idOuterOptimMetsConcRates;
    case {'init', 'local'}
        nOptimMetsInitConcs = nNonPoolMets;
        idOptimMetsInitConcs =idNonPoolMets;
        nOptimMetsConcRates = nNonPoolMets;
        idOptimMetsConcRates =idNonPoolMets;
end


%% Bounds for initConcs, switchTimeConcRates, switchTimeFluxes, and switch times
lb.initConcs= zeros(nOptimMetsInitConcs, 1)+optionsMFA.lb.initConcs;
lb.switchTimeConcRates= zeros(nOptimMetsConcRates, nSwitchTimes+2)+optionsMFA.lb.switchTimeConcRates;
lb.switchTimeFluxes = zeros(nIndFluxes, nSwitchTimes+2)+optionsMFA.lb.switchTimeFluxes;
lb.switchTimes =  zeros(nSwitchTimes,1)+optionsMFA.lb.switchTimes+optionsMFA.minKnotTimeDiff;

ub.initConcs= zeros(nOptimMetsInitConcs, 1)+optionsMFA.ub.initConcs;
ub.switchTimeConcRates= zeros(nOptimMetsConcRates, nSwitchTimes+2)+optionsMFA.ub.switchTimeConcRates;
ub.switchTimeFluxes = zeros(nIndFluxes, nSwitchTimes+2)+optionsMFA.ub.switchTimeFluxes;
ub.switchTimes =  zeros(nSwitchTimes,1)+optionsMFA.ub.switchTimes-optionsMFA.minKnotTimeDiff;

switch optimType
    case {'metaheuristic'}
        if optionsMFA.isIncludeAllInitConcsInParamMH
            loc = isEvalConcMets(1:nNonPoolMets);
            lb.initConcs(loc,:) = optionsMFA.lb.mesInitConcs;
            ub.initConcs(loc,:) = optionsMFA.ub.mesInitConcs;
        end
    case {'init', 'local'}
        loc = isEvalConcMets(1:nNonPoolMets);
        lb.initConcs(loc,:) = optionsMFA.lb.mesInitConcs;
        lb.switchTimeConcRates(loc,:) = optionsMFA.lb.mesKnotConcRates;
        ub.initConcs(loc,:) = optionsMFA.ub.mesInitConcs;
        ub.switchTimeConcRates(loc,:) = optionsMFA.ub.mesKnotConcRates;
end


%% Expand parameter bound according to expandRatio
lb.initConcs = lb.initConcs/expandRatio;
lb.switchTimeConcRates= lb.switchTimeConcRates*expandRatio;
lb.switchTimeFluxes = lb.switchTimeFluxes/expandRatio;
lb.switchTimes =  lb.switchTimes -optionsMFA.minKnotTimeDiff+optionsMFA.minKnotTimeDiff/expandRatio;

ub.initConcs= ub.initConcs*expandRatio;
ub.switchTimeConcRates= ub.switchTimeConcRates*expandRatio;
ub.switchTimeFluxes = ub.switchTimeFluxes*expandRatio;
ub.switchTimes =  ub.switchTimes +optionsMFA.minKnotTimeDiff-optionsMFA.minKnotTimeDiff/expandRatio;

%% Bounds of fluxes
switch optimType
    case {'init','local'}
        lb.switchTimeFluxesAll = zeros(nRxns,nSwitchTimes+2)+optionsMFA.lb.switchTimeFluxes;
        ub.switchTimeFluxesAll = zeros(nRxns,nSwitchTimes+2)+optionsMFA.ub.switchTimeFluxes;
        lb.switchTimeFluxesAll = lb.switchTimeFluxesAll/expandRatio;
        ub.switchTimeFluxesAll = ub.switchTimeFluxesAll*expandRatio;
end


end
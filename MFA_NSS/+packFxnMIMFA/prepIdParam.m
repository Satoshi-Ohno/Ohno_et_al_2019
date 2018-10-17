%% Preparation of index of parameters
function idParam = prepIdParam(model, expData, optionsMFA, optimType)

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
        lb = optionsMFA.lbMH;
        ub = optionsMFA.ubMH;
    case {'init', 'local'}
        nOptimMetsInitConcs = nNonPoolMets;
        idOptimMetsInitConcs =idNonPoolMets;
        nOptimMetsConcRates = nNonPoolMets;
        idOptimMetsConcRates =idNonPoolMets;
        lb = optionsMFA.lbLocal;
        ub = optionsMFA.ubLocal;
end


%% ID of parameters
nParam = 0;

idParam.initConcs = nParam+(1:nOptimMetsInitConcs);
nParam = nParam + length(idParam.initConcs);

idParam.switchTimeConcRates = nParam+(1:nOptimMetsConcRates*(nSwitchTimes+2));
nParam = nParam + length(idParam.switchTimeConcRates);

idParam.switchTimeFluxes = nParam+(1:nIndFluxes * (nSwitchTimes+2));
nParam = nParam + length(idParam.switchTimeFluxes);

idParam.switchTimes = nParam + (1:nSwitchTimes);
nParam = nParam + length(idParam.switchTimes);

idParam.nParam =nParam;
end

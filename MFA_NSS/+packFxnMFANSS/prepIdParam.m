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
        nOptimMetsConcs =  nMHMetsConcs;
        idOptimMetsConcs =idMHMetsConcs;
        nOptimMetsConcRates = nMHMetsConcRates;
        idOptimMetsConcRates = idMHMetsConcRates;
        lb = optionsMFA.lbMH;
        ub = optionsMFA.ubMH;
    case {'init', 'local'}
        nOptimMetsConcs = nNonPoolMets;
        idOptimMetsConcs =idNonPoolMets;
        nOptimMetsConcRates = nNonPoolMets;
        idOptimMetsConcRates =idNonPoolMets;
        lb = optionsMFA.lbLocal;
        ub = optionsMFA.ubLocal;
end


%% ID of parameters
nParam = 0;

if optionsMFA.isUseConcAsParam
    idParam.initConcs = nParam+(1:nOptimMetsConcs);
    idParam.concs = nParam+(1:nOptimMetsConcs*(nSwitchTimes+2));
    nParam = nParam + length(idParam.concs);
    
    idParam.initConcRates = nParam+(1:nOptimMetsConcRates);
    idParam.concRates = nParam+(1:nOptimMetsConcRates);
    nParam = nParam + length(idParam.concRates);
else
    idParam.initConcs = nParam+(1:nOptimMetsConcs);
    idParam.concs = nParam+(1:nOptimMetsConcs);
    nParam = nParam + length(idParam.concs);
    
    idParam.initConcRates = nParam+(1:nOptimMetsConcRates);
    idParam.concRates= nParam+(1:nOptimMetsConcRates*(nSwitchTimes+2));
    nParam = nParam + length(idParam.concRates);
end

idParam.fluxes = nParam+(1:nIndFluxes * (nSwitchTimes+2));
nParam = nParam + length(idParam.fluxes);

idParam.switchTimes = nParam + (1:nSwitchTimes);
nParam = nParam + length(idParam.switchTimes);

idParam.nParam =nParam;
end

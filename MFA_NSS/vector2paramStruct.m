function    paramStruct = vector2paramStruct(param, model, expData, optionsMFA)

nRxns = optionsMFA.varSet.nRxns;
nNonPoolMets = optionsMFA.varSet.nNonPoolMets;
nSwitchTimes = optionsMFA.varSet.nSwitchTimes;

idParam =optionsMFA.idParamLocal;
transformMat = optionsMFA.transformMat;

%% switch time
switchTimes = param(idParam.switchTimes);
switchTimes = sort(switchTimes);

%% flux at switch time
switchTimeFluxes = transformMat.param2Flux * param;
switchTimeFluxes = reshape(switchTimeFluxes, nRxns, nSwitchTimes+2);

%% initial metabolite concentrations
initConcs = transformMat.param2InitConc * param;

%% metabolite concentrations at switch time
switchTimeConcs = transformMat.param2Conc * param;
switchTimeConcs = reshape(switchTimeConcs, nNonPoolMets, nSwitchTimes+2);

%% metabolite concentrations at switch time
initConcRates = transformMat.param2InitConcRate * param;

%% change in metabolite concentrations at switch time
switchTimeConcRates = transformMat.param2ConcRate * param;
switchTimeConcRates = reshape(switchTimeConcRates, nNonPoolMets, nSwitchTimes+2);

%% polynomial coefficient of metabolite concetnrations
pCoefConcs = transformMat.param2pCoefConc*param;
pCoefConcs = reshape(pCoefConcs, nNonPoolMets, 3*(nSwitchTimes+1));

%% polynomial coefficient of flux
pCoefFluxes = transformMat.param2pCoefFlux * param;
pCoefFluxes = reshape(pCoefFluxes,nRxns, 2*(nSwitchTimes+1));

%% merge
paramStruct.initConcs = initConcs;
paramStruct.switchTimeConcs = switchTimeConcs;
paramStruct.initConcRates = initConcRates;
paramStruct.switchTimeConcRates = switchTimeConcRates;
paramStruct.switchTimeFluxes = switchTimeFluxes;
paramStruct.switchTimes= switchTimes;
paramStruct.pCoefConcs= pCoefConcs;
paramStruct.pCoefFluxes= pCoefFluxes;


end

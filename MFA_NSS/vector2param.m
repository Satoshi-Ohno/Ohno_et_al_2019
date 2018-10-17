function    [switchTimeFluxes, initConcs, switchTimeConcRates, switchTimes] = ...
    vector2param(param, model, expData, optionsMFA)


field2var(optionsMFA.varSet)

        idParam =optionsMFA.idParamLocal;
convertMat = optionsMFA.convertMat;

%% switch time
switchTimes = param(idParam.switchTimes);
switchTimes = sort(switchTimes);

%% flux at switch time
switchTimeFluxes = convertMat.param2Flux * param;
switchTimeFluxes = reshape(switchTimeFluxes, nRxns, nSwitchTimes+2);

%% initial metabolite concentrations
initConcs = convertMat.param2InitConc * param;

%% change in metabolite concentrations at switch time
switchTimeConcRates = convertMat.param2ConcRate * param;
switchTimeConcRates = reshape(switchTimeConcRates, nNonPoolMets, nSwitchTimes+2);

end

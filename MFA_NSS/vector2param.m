function    [knotFluxes, initConcs, knotConcRates, knots, coefCorrCompt] = ...
    vector2param170807(param, model, expData, optionsMFA)


field2var(optionsMFA.varSet)

        idParam =optionsMFA.idParamLocal;
convertMat = optionsMFA.convertMat;

%% knot time
knots = param(idParam.knots);
knots = sort(knots);
% fullKnots = [0, knots', expData.time(end)];

%% coefCorrCompt
coefCorrCompt = param(idParam.coefCorrCompt);

%% flux at knot time
knotFluxes = convertMat.param2Flux * param;
knotFluxes = reshape(knotFluxes, nRxns, nKnots+2);

%% initial metabolite concentrations
initConcs = convertMat.param2InitConc * param;

%% change in metabolite concentrations at knot time
knotConcRates = convertMat.param2ConcRate * param;
knotConcRates = reshape(knotConcRates, nNonPoolMets, nKnots+2);

end

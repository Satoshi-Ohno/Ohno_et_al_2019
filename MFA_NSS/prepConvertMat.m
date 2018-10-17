%% Make matrix to obtain variales from parameter vector
function convertMat = prepConvertMat170803(model, expData, optionsMFA)

field2var(optionsMFA.varSet)
idParamLocal = optionsMFA.idParamLocal;
nParamLocal = idParamLocal.nParam;
idParamMH = optionsMFA.idParamMH;
nParamMH = idParamMH.nParam;


%% param -> initConc
matParam2InitConc  = sparse(nNonPoolMets, nParamLocal);
row = 1:nNonPoolMets;
col =  optionsMFA.idParamLocal.initConcs;
matParam2InitConc(row, col) = matParam2InitConc(row, col) + eye(nNonPoolMets);

%% param -> concRate
matParam2ConcRate = sparse(nNonPoolMets*(nSwitchTimes+2), nParamLocal);
for k = 1 : nSwitchTimes+2
    row = nNonPoolMets*(k-1)+(1:nNonPoolMets);
    col =  idParamLocal.switchTimeConcRates(nNonPoolMets*(k-1)+(1:nNonPoolMets));
    matParam2ConcRate(row, col) = matParam2ConcRate(row, col) + eye(nNonPoolMets);
end

%% param -> indFlux
matParam2IndFlux = sparse(nIndFluxes*(nSwitchTimes+2), nParamLocal);
for k = 1 : nSwitchTimes+2
    row = nIndFluxes*(k-1)+(1:nIndFluxes);
    col =  idParamLocal.switchTimeFluxes(nIndFluxes*(k-1)+(1:nIndFluxes));
    matParam2IndFlux(row, col) = matParam2IndFlux(row, col) + eye(nIndFluxes);
end

%% param -> concRate + indFlux
matParam2ConcRateIndFlux = sparse((nNonPoolMets+nIndFluxes)*(nSwitchTimes+2),nParamLocal);
for k = 1 : nSwitchTimes+2
    row0 = (nNonPoolMets+nIndFluxes)*(k-1) + (1:(nNonPoolMets+nIndFluxes));
    row1 = (nNonPoolMets)*(k-1) + (1:(nNonPoolMets));
    row2 = (nIndFluxes)*(k-1) + (1:(nIndFluxes));
    matParam2ConcRateIndFlux(row0,:) = [matParam2ConcRate(row1,:);matParam2IndFlux(row2,:)];
end

%% concRate + indFlux -> flux
matConcRateIndFlux2Flux = sparse(nRxns*(nSwitchTimes+2), (nNonPoolMets+nIndFluxes)*(nSwitchTimes+2));
tmpMat = full(model.invS);
for k = 1 : nSwitchTimes+2
    row = nRxns*(k-1) + (1:nRxns);
    col = (nNonPoolMets+nIndFluxes)*(k-1) + (1:(nNonPoolMets+nIndFluxes));
    matConcRateIndFlux2Flux(row, col) = matConcRateIndFlux2Flux(row, col) + tmpMat;
end

%% flux -> netFlux
matFlux2NetFlux = sparse(nNetRxns*(nSwitchTimes+2), nRxns*(nSwitchTimes+2));
tmpMat = model.matRaw2Net;
for k = 1 : nSwitchTimes+2
    row = nNetRxns*(k-1)+(1:nNetRxns);
    col =  nRxns*(k-1)+(1:nRxns);
    matFlux2NetFlux(row, col) = matFlux2NetFlux(row, col)+tmpMat;
end

%% flux -> concRate
matFlux2ConcRate= sparse(nNonPoolMets*(nSwitchTimes+2), nRxns*(nSwitchTimes+2));
tmpMat = full(model.S(1:nNonPoolMets,:));
for k = 1 : nSwitchTimes+2
    row = nNonPoolMets*(k-1)+(1:nNonPoolMets);
    col =  nRxns*(k-1)+(1:nRxns);
    matFlux2ConcRate(row, col) = matFlux2ConcRate(row, col)+tmpMat;
end

%% flux -> concRateAllMets
matFlux2ConcRateAllMet = sparse(nMets*(nSwitchTimes+2), nRxns*(nSwitchTimes+2));
tmpMat = full(model.S);
for k = 1 : nSwitchTimes+2
    row = nMets*(k-1)+(1:nMets);
    col =  nRxns*(k-1)+(1:nRxns);
    matFlux2ConcRateAllMet(row, col) = matFlux2ConcRateAllMet(row, col)+tmpMat;
end

%% param -> flux, netFlux, concRates, concRatesAllMet
matParam2Flux = matConcRateIndFlux2Flux*matParam2ConcRateIndFlux;
matParam2NetFlux = matFlux2NetFlux*matConcRateIndFlux2Flux*matParam2ConcRateIndFlux;
matParam2ConcRateAllMet = matFlux2ConcRateAllMet*matConcRateIndFlux2Flux*matParam2ConcRateIndFlux;

%% paramMH -> paramLocal 
if optionsMFA.isUseQPInMH
    matParamMH2ParamLocal = sparse(nParamLocal, nParamMH);
    tmpFieldNames = fieldnames(idParamMH);
    for f = 1 : length(tmpFieldNames)
        switch tmpFieldNames{f}
            case {'switchTimeFluxes', 'switchTimes'}
                row = idParamLocal.(tmpFieldNames{f});
                col = idParamMH.(tmpFieldNames{f});
            case {'initConcs'}
                row = idParamLocal.(tmpFieldNames{f})(idOuterOptimMetsInitConcs);
                col = idParamMH.(tmpFieldNames{f});
            case {'switchTimeConcRates'}
                isConcRateOuterOptimMets = false(nNonPoolMets, nSwitchTimes+2);
                isConcRateOuterOptimMets(idOuterOptimMetsConcRates,:) = true;
                isConcRateOuterOptimMets = reshape(isConcRateOuterOptimMets, nNonPoolMets*(nSwitchTimes+2), 1);
                row = idParamLocal.switchTimeConcRates(isConcRateOuterOptimMets);
                col = idParamMH.(tmpFieldNames{f});
            otherwise
                row = [];
                col = [];
        end
        matParamMH2ParamLocal(row, col) = matParamMH2ParamLocal(row, col) + eye(length(row));
    end
else 
    matParamMH2ParamLocal = speye(nParamLocal, nParamMH);
end

%% Identify parameters in log scale
isLogParamMH = false(nParamMH,1);
idLog = [idParamMH.initConcs, idParamMH.switchTimeFluxes, idParamMH.switchTimes];
isLogParamMH(idLog) = true;

%% paramLocal -> paramMH
matParamLocal2ParamMH = matParamMH2ParamLocal';

%% paramMH -> paramMHInd
if isfield(optionsMFA, 'isIndParams') && isfield(optionsMFA.isIndParams, 'MH')
    idIndParamsMH = find(optionsMFA.isIndParams.MH);
    matParamMH2Ind = sparse(length(idIndParamsMH), nParamMH);
    matParamMH2Ind(:, idIndParamsMH) = speye(length(idIndParamsMH));
else
    matParamMH2Ind = [];
end

%% paramMH -> concRate
if optionsMFA.isUseQPInMH
    matParamMH2ConcRate = sparse(nOuterOptimMetsConcRates*(nSwitchTimes+2), nParamMH);
    for k = 1 : nSwitchTimes+2
        row = nOuterOptimMetsConcRates*(k-1)+(1:nOuterOptimMetsConcRates);
        col =  idParamMH.switchTimeConcRates(nOuterOptimMetsConcRates*(k-1)+(1:nOuterOptimMetsConcRates));
        matParamMH2ConcRate(row, col) = matParamMH2ConcRate(row, col) + eye(nOuterOptimMetsConcRates);
    end
else 
    matParamMH2ConcRate = matParam2ConcRate;
end

%% paramMH -> indFlux
    matParamMH2IndFlux = sparse(nIndFluxes*(nSwitchTimes+2), nParamMH);
    for k = 1 : nSwitchTimes+2
        row = nIndFluxes*(k-1)+(1:nIndFluxes);
        col =  idParamMH.switchTimeFluxes(nIndFluxes*(k-1)+(1:nIndFluxes));
        matParamMH2IndFlux(row, col) = matParamMH2IndFlux(row, col) + eye(nIndFluxes);
    end
    
%% flux -> nonIndFlux
matFlux2NonIndFlux = sparse(nNonIndFluxes*(nSwitchTimes+2), nRxns*(nSwitchTimes+2));
for k = 1 : nSwitchTimes+2
    row = nNonIndFluxes*(k-1)+(1:nNonIndFluxes);
    col =  nRxns*(k-1)+idNonIndFluxes;
    matFlux2NonIndFlux(row, col) = matFlux2NonIndFlux(row, col) + eye(nNonIndFluxes);
end

%% param -> nonIndFlux
matParam2NonIndFlux = matFlux2NonIndFlux*matParam2Flux;

%% constraint for mass balance

matCalcParam2ConcRate = matFlux2ConcRate*matConcRateIndFlux2Flux*matParam2ConcRateIndFlux;
matConstrMassBalance = matParam2ConcRate - matCalcParam2ConcRate;
matConstrMassBalance(abs(matConstrMassBalance)<=10^-9) =0;
isRedundantConstr = all(matConstrMassBalance==0,2);
matConstrMassBalance = matConstrMassBalance(~isRedundantConstr, :);

%% merge
convertMat.param2InitConc = matParam2InitConc;
convertMat.param2IndFlux = matParam2IndFlux;
convertMat.param2ConcRateIndFlux = matParam2ConcRateIndFlux;
convertMat.param2Flux = matParam2Flux;
convertMat.param2NonIndFlux = matParam2NonIndFlux;
convertMat.param2NetFlux = matParam2NetFlux;
convertMat.param2ConcRate = matParam2ConcRate;
convertMat.param2ConcRateAllMet = matParam2ConcRateAllMet;
convertMat.paramMH2ParamLocal = matParamMH2ParamLocal;
convertMat.paramLocal2ParamMH = matParamLocal2ParamMH;
convertMat.paramMH2Ind = matParamMH2Ind;
convertMat.paramMH2ConcRate = matParamMH2ConcRate;
convertMat.paramMH2IndFlux = matParamMH2IndFlux;
convertMat.constrMassBalance = matConstrMassBalance;

tmpFieldNames = fieldnames(convertMat);
for f = 1 : length(tmpFieldNames)
    tmpMat = convertMat.(tmpFieldNames{f});
    tmpMat(abs(tmpMat)<=10^-9) = 0;  % remove error
    convertMat.(tmpFieldNames{f})= tmpMat;
end

convertMat.isLogParamMH = isLogParamMH;

end




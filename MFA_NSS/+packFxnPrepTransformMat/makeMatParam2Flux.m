%% Make transformation matrix (param -> flux)
function matParam2Flux = ...
    makeMatParam2Flux(model, expData, optionsMFA, fullSwitchTimes, inputFlag)

%% prepare variables
% idParamLocal = optionsMFA.idParamLocal;
% nParamLocal = idParamLocal.nParam;
% idParamMH = optionsMFA.idParamMH;
% nParamMH = idParamMH.nParam;
% nMets = optionsMFA.varSet.nMets;
% nRxns = optionsMFA.varSet.nRxns;
% nNonPoolMets = optionsMFA.varSet.nNonPoolMets;
% idNonPoolMets = optionsMFA.varSet.idNonPoolMets;
% nSwitchTimes = optionsMFA.varSet.nSwitchTimes;
% idEvalConcMets = optionsMFA.varSet.idEvalConcMets;
% nIndFluxes = optionsMFA.varSet.nIndFluxes;
% idNonIndFluxes = optionsMFA.varSet.idNonIndFluxes;
% nNonIndFluxes = optionsMFA.varSet.nNonIndFluxes;
% nNetRxns = optionsMFA.varSet.nNetRxns;
% nMHMetsConcRates = optionsMFA.varSet.nMHMetsConcRates;
% idMHMetsConcs = optionsMFA.varSet.idMHMetsConcs;
% idMHMetsConcRates = optionsMFA.varSet.idMHMetsConcRates;

%% make matrix
% inputFlag: 
% 11, isempty(fullSwitchTimes) == true & optionsMFA.isUseConcAsParam == true
% 12, isempty(fullSwitchTimes) == true & optionsMFA.isUseConcAsParam == false
% 21, isempty(fullSwitchTimes) == false & optionsMFA.isUseConcAsParam == true
% 22, isempty(fullSwitchTimes) == false & optionsMFA.isUseConcAsParam == false

switch inputFlag
    case 11
        matParam2Flux = [];
    case {12, 21}
        matParam2ConcRateIndFlux = ...
            makeMatParam2ConcRateIndFlux(model, expData, optionsMFA, fullSwitchTimes, inputFlag);
        matConcRateIndFlux2Flux = ...
            makeMatConcRateIndFlux2Flux(model, expData, optionsMFA, fullSwitchTimes, inputFlag);            
        matParam2Flux = matConcRateIndFlux2Flux*matParam2ConcRateIndFlux;
    case 22
        matParam2Flux = optionsMFA.transformMat.param2Flux;
end


end

%% Make convert matrix (param -> concRate + indFlux)
function matParam2ConcRateIndFlux = ...
    makeMatParam2ConcRateIndFlux(model, expData, optionsMFA, fullSwitchTimes, inputFlag)

%% prepare variables
idParamLocal = optionsMFA.idParamLocal;
nParamLocal = idParamLocal.nParam;
% idParamMH = optionsMFA.idParamMH;
% nParamMH = idParamMH.nParam;
% nMets = optionsMFA.varSet.nMets;
% nRxns = optionsMFA.varSet.nRxns;
nNonPoolMets = optionsMFA.varSet.nNonPoolMets;
% idNonPoolMets = optionsMFA.varSet.idNonPoolMets;
nSwitchTimes = optionsMFA.varSet.nSwitchTimes;
% idEvalConcMets = optionsMFA.varSet.idEvalConcMets;
nIndFluxes = optionsMFA.varSet.nIndFluxes;
% idNonIndFluxes = optionsMFA.varSet.idNonIndFluxes;
% nNonIndFluxes = optionsMFA.varSet.nNonIndFluxes;
% nNetRxns = optionsMFA.varSet.nNetRxns;
% nMHMetsConcRates = optionsMFA.varSet.nMHMetsConcRates;
% idMHMetsConcs = optionsMFA.varSet.idMHMetsConcs;
% idMHMetsConcRates = optionsMFA.varSet.idMHMetsConcRates;

%% make matrix

matParam2ConcRate = optionsMFA.matParam2ConcRate;
matParam2IndFlux = optionsMFA.matParam2IndFlux;
matParam2ConcRateIndFlux = sparse((nNonPoolMets+nIndFluxes)*(nSwitchTimes+2),nParamLocal);
for k = 1 : nSwitchTimes+2
    row0 = (nNonPoolMets+nIndFluxes)*(k-1) + (1:(nNonPoolMets+nIndFluxes));
    row1 = (nNonPoolMets)*(k-1) + (1:(nNonPoolMets));
    row2 = (nIndFluxes)*(k-1) + (1:(nIndFluxes));
    matParam2ConcRateIndFlux(row0,:) = [matParam2ConcRate(row1,:);matParam2IndFlux(row2,:)];
end


end

%% Make convert matrix (concRate + indFlux -> flux)
function matConcRateIndFlux2Flux = ...
    makeMatConcRateIndFlux2Flux(model, expData, optionsMFA, fullSwitchTimes, inputFlag)

%% prepare variables
% idParamLocal = optionsMFA.idParamLocal;
% nParamLocal = idParamLocal.nParam;
% idParamMH = optionsMFA.idParamMH;
% nParamMH = idParamMH.nParam;
% nMets = optionsMFA.varSet.nMets;
nRxns = optionsMFA.varSet.nRxns;
nNonPoolMets = optionsMFA.varSet.nNonPoolMets;
% idNonPoolMets = optionsMFA.varSet.idNonPoolMets;
nSwitchTimes = optionsMFA.varSet.nSwitchTimes;
% idEvalConcMets = optionsMFA.varSet.idEvalConcMets;
nIndFluxes = optionsMFA.varSet.nIndFluxes;
% idNonIndFluxes = optionsMFA.varSet.idNonIndFluxes;
% nNonIndFluxes = optionsMFA.varSet.nNonIndFluxes;
% nNetRxns = optionsMFA.varSet.nNetRxns;
% nMHMetsConcRates = optionsMFA.varSet.nMHMetsConcRates;
% idMHMetsConcs = optionsMFA.varSet.idMHMetsConcs;
% idMHMetsConcRates = optionsMFA.varSet.idMHMetsConcRates;

%% make matrix

matConcRateIndFlux2Flux = sparse(nRxns*(nSwitchTimes+2), (nNonPoolMets+nIndFluxes)*(nSwitchTimes+2));
tmpMat = full(model.invS);
for k = 1 : nSwitchTimes+2
    row = nRxns*(k-1) + (1:nRxns);
    col = (nNonPoolMets+nIndFluxes)*(k-1) + (1:(nNonPoolMets+nIndFluxes));
    matConcRateIndFlux2Flux(row, col) = matConcRateIndFlux2Flux(row, col) + tmpMat;
end


end
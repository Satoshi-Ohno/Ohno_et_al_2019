%% Make transformation matrix (param -> concExpTime)
function matParam2pCoefFlux = ...
    makeMatParam2pCoefFlux(model, expData, optionsMFA, fullSwitchTimes, inputFlag)

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
    case {11, 12}
        matParam2pCoefFlux = [];
    case {21,22}
        matFlux2pCoefFlux = ...
            makeMatFlux2pCoefFlux(model, expData, optionsMFA, fullSwitchTimes);        
        matParam2pCoefFlux = matFlux2pCoefFlux * optionsMFA.matParam2Flux;
end

end

%% concConst + concRate -> pCoefFlux
function matFlux2pCoefFlux = ...
    makeMatFlux2pCoefFlux(model, expData, optionsMFA, fullSwitchTimes)

nRxns = optionsMFA.varSet.nRxns;
nSwitchTimes = optionsMFA.varSet.nSwitchTimes;

matFlux2pCoefFlux = sparse(...
    nRxns*2*(nSwitchTimes+1),...
    nRxns*(nSwitchTimes+2)...
    );
for k = 1 : nSwitchTimes+1
    tmpMatCoef = zeros(2,2);
    tmpMatCoef(1,1) = fullSwitchTimes(k+1);
    tmpMatCoef(1,2) = -fullSwitchTimes(k);
    tmpMatCoef(2,1) = -1;
    tmpMatCoef(2,2) = 1;
    
    tmpMatIntegral = eye(2);
    
    tmpMat = zeros(size(tmpMatIntegral,1)*nRxns, size(tmpMatIntegral,2)*nRxns);
    for r = 1 : nRxns
        row = size(tmpMatIntegral,1)*(r-1)+(1:size(tmpMatIntegral,1));
        col = size(tmpMatIntegral,2)*(r-1)+(1:size(tmpMatIntegral,2));
        tmpMat(row, col) = tmpMat(row, col) + tmpMatIntegral * tmpMatCoef;
    end
    
    row = nRxns*2*(k-1) ...
        + repmat((1:nRxns)', 1, 2) ...
        + repmat(nRxns*(0:1), nRxns,  1);
    row = reshape(row', numel(row), 1);
    col = repmat((1:nRxns)', 1, 2) + repmat(nRxns*(k-1:k), nRxns, 1);
    col = reshape(col', numel(col), 1);
    
    matFlux2pCoefFlux(row, col) = matFlux2pCoefFlux(row, col) + ...
        1/(fullSwitchTimes(k+1)-fullSwitchTimes(k)) * sparse(tmpMat);
end


end


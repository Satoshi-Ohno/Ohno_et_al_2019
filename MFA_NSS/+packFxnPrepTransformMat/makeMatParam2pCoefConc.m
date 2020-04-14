%% Make transformation matrix (param -> concExpTime)
function matParam2pCoefConc = ...
    makeMatParam2pCoefConc(model, expData, optionsMFA, fullSwitchTimes, inputFlag)

%% prepare variables
% idParamLocal = optionsMFA.idParamLocal;
% nParamLocal = idParamLocal.nParam;
% idParamMH = optionsMFA.idParamMH;
% nParamMH = idParamMH.nParam;
% nMets = optionsMFA.varSet.nMets;
% nRxns = optionsMFA.varSet.nRxns;
% nNonPoolMets = optionsMFA.varSet.nNonPoolMets;
idNonPoolMets = optionsMFA.varSet.idNonPoolMets;
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
        matParam2pCoefConc = [];
    case {21,22}
        matParam2ConcConst = optionsMFA.matParam2ConcConst;
        matParam2ConcRate = optionsMFA.matParam2ConcRate;
        matConcConstRate2pCoefConc = ...
            makeMatConcConstRate2pCoefConc(model, expData, optionsMFA, fullSwitchTimes, idNonPoolMets);
        matParam2pCoefConc = matConcConstRate2pCoefConc * [matParam2ConcConst;matParam2ConcRate];
end

end

%% concConst + concRate -> pCoefConc
function matConcConstRate2pCoefConc = ...
    makeMatConcConstRate2pCoefConc(model, expData, optionsMFA, fullSwitchTimes, idMetsAll)

nNonPoolMets = optionsMFA.varSet.nNonPoolMets;
nSwitchTimes = optionsMFA.varSet.nSwitchTimes;
idMetsAll = colVec(idMetsAll);

matConcConstRate2pCoefConc = sparse(...
    length(idMetsAll)*3*(nSwitchTimes+1),...
    nNonPoolMets*((nSwitchTimes+1) + (nSwitchTimes+2))...
    );
for k = 1 : nSwitchTimes+1
    tmpMatCoef = zeros(3,3);
    tmpMatCoef(1,1) = fullSwitchTimes(k+1)-fullSwitchTimes(k);
    tmpMatCoef(2,2) = fullSwitchTimes(k+1);
    tmpMatCoef(2,3) = -fullSwitchTimes(k);
    tmpMatCoef(3,2) = -1;
    tmpMatCoef(3,3) = 1;
    
    tmpMatIntegral = eye(3);
    tmpMatIntegral(3,3) = 1/2;
    
    tmpMat = zeros(size(tmpMatIntegral,1)*length(idMetsAll), size(tmpMatIntegral,2)*length(idMetsAll));
    for m = 1 : length(idMetsAll)
        row = size(tmpMatIntegral,1)*(m-1) + (1:size(tmpMatIntegral,1));
        col = size(tmpMatIntegral,2)*(m-1) + (1:size(tmpMatIntegral,2));
        tmpMat(row, col) = tmpMat(row, col) + tmpMatIntegral * tmpMatCoef;
    end
    
    
    row = length(idMetsAll)*3*(k-1) ...
        + repmat(idMetsAll,1, 3) ...
        + repmat(length(idMetsAll)*(0:2),length(idMetsAll), 1);
    row = reshape(row', numel(row), 1);
    
    colConcConst = nNonPoolMets*(k-1)+idMetsAll;
    colConcRate1 = nNonPoolMets*(nSwitchTimes+1)+nNonPoolMets*(k-1)+idMetsAll;
    colConcRate2 = nNonPoolMets*(nSwitchTimes+1)+nNonPoolMets*(k)+idMetsAll;
    col = [colConcConst, colConcRate1, colConcRate2];
    col = reshape(col', numel(col), 1);
    matConcConstRate2pCoefConc(row, col) = matConcConstRate2pCoefConc(row, col) + ...
        1/(fullSwitchTimes(k+1)-fullSwitchTimes(k)) * sparse(tmpMat);
    
end

end



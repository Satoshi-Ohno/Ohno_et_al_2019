%% Make matrix to obtain variales from parameter vector and switch time
function convertMat = prepConvertMatDepSwitchTimes(model, expData, optionsMFA, fullSwitchTimes)
idNonPoolMets = optionsMFA.varSet.idNonPoolMets;
idEvalConcMets = optionsMFA.varSet.idEvalConcMets;
nSwitchTimes = optionsMFA.varSet.nSwitchTimes;

convertMat = optionsMFA.convertMat;
idParam = optionsMFA.idParamLocal;
nParam = idParam.nParam;
timeConstrQP = linspace(0, expData.time(end), optionsMFA.nStepTimeConstrQP);

idTimeInterval = zeros(1,length(timeConstrQP));
for k = 1 : nSwitchTimes+1
    loc = timeConstrQP>=fullSwitchTimes(k) & timeConstrQP<=fullSwitchTimes(k+1);
    idTimeInterval(loc) = k;    
end

%% param -> concConst
matParam2ConcConst =  ...
    makeMatParam2ConcConst(model, expData, optionsMFA, fullSwitchTimes);

%% param -> concRate
matParam2ConcRate=convertMat.param2ConcRate;

%% concConst + concRates -> concs
matConcConstRate2Conc = ...
    makeMatConcConstRate2Conc(model, expData, optionsMFA, fullSwitchTimes, timeConstrQP, idNonPoolMets);

%% concConst + concRates -> concsExpTime
matConcConstRate2ConcExpTime = ...
    makeMatConcConstRate2Conc(model, expData, optionsMFA, fullSwitchTimes, expData.time, idEvalConcMets);

%% paramLocal -> conc

matParam2Conc = matConcConstRate2Conc * [matParam2ConcConst;matParam2ConcRate];
matParam2ConcExpTime = matConcConstRate2ConcExpTime * [matParam2ConcConst;matParam2ConcRate];


%% param -> pCoefConcs 

matConcConstRate2pCoefConc = ...
    makeMatConcConstRate2pCoefConc(model, expData, optionsMFA, fullSwitchTimes, idNonPoolMets);

matParam2pCoefConc = matConcConstRate2pCoefConc * [matParam2ConcConst;matParam2ConcRate];

%% paramLocal -> pCoefFluxes
matFlux2pCoefFlux = ...
    makeMatFlux2pCoefFlux(model, expData, optionsMFA, fullSwitchTimes);

matParam2pCoefFlux = matFlux2pCoefFlux * convertMat.param2Flux;

%% Constraint for reaction dependent number of time intervvals
varType = 'concRate';
matConcRate2NotForOptim = ...
    makeVar2NotForOptim(model, expData, optionsMFA, fullSwitchTimes, varType);
matCalcConcRate2NotForOptim = ...
    makeCalcVar2NotForOptim(model, expData, optionsMFA, fullSwitchTimes, varType);

varType = 'flux';
matFlux2NotForOptim = ...
    makeVar2NotForOptim(model, expData, optionsMFA, fullSwitchTimes, varType);
matCalcFlux2NotForOptim = ...
    makeCalcVar2NotForOptim(model, expData, optionsMFA, fullSwitchTimes, varType);


matConstrConcRateNotForOptim = ...
    (matCalcConcRate2NotForOptim - matConcRate2NotForOptim) * convertMat.param2ConcRate;
matConstrFluxNotForOptim = ...
    (matCalcFlux2NotForOptim - matFlux2NotForOptim) * convertMat.param2Flux;


%% merge
preConvertMat.param2ConcConst = matParam2ConcConst;
preConvertMat.param2Conc = matParam2Conc;
preConvertMat.param2ConcExpTime = matParam2ConcExpTime;
% preConvertMat.concExpTime2ConcMerged = matConcExpTime2ConcMerged;
preConvertMat.param2pCoefConc = matParam2pCoefConc;
preConvertMat.param2pCoefFlux = matParam2pCoefFlux;
preConvertMat.constrConcRateNotForOptim = matConstrConcRateNotForOptim;
preConvertMat.constrFluxNotForOptim = matConstrFluxNotForOptim;

tmpFieldNames = fieldnames(preConvertMat);
for f = 1 : length(tmpFieldNames)
    tmpMat = preConvertMat.(tmpFieldNames{f});
    convertMat.(tmpFieldNames{f})= tmpMat;
end

end

%% concConst + concRate -> conc
function matConcConstRate2Conc = ...
    makeMatConcConstRate2Conc(model, expData, optionsMFA, fullSwitchTimes, timeConstr, idMetsAll, isEarlierTI)

nNonPoolMets = optionsMFA.varSet.nNonPoolMets;
nSwitchTimes = optionsMFA.varSet.nSwitchTimes;

idTimeInterval = zeros(1,length(timeConstr));
    switchTimeOrder = 1 : nSwitchTimes+1;
for k = switchTimeOrder
    loc = timeConstr>=fullSwitchTimes(k) & timeConstr<=fullSwitchTimes(k+1);
    idTimeInterval(loc) = k;
end


matConcConstRate2Conc = sparse(...
    length(timeConstr)*length(idMetsAll),...
    nNonPoolMets*((nSwitchTimes+1) + (nSwitchTimes+2))...
    ); 
kk = 0;
for k = 1 : nSwitchTimes+1 
    tmpMatCoef = zeros(3,3);
    tmpMatCoef(1,1) = fullSwitchTimes(k+1)-fullSwitchTimes(k);
    tmpMatCoef(2,2) = fullSwitchTimes(k+1);
    tmpMatCoef(2,3) = -fullSwitchTimes(k);
    tmpMatCoef(3,2) = -1;
    tmpMatCoef(3,3) = 1;
    
    loc = idTimeInterval == k;
    tmpTime = timeConstr(loc);
    if isempty(tmpTime)
        continue
    end
    tmpMatTime = [ones(length(tmpTime), 1),  tmpTime',  1/2*tmpTime'.^2 ];
    
    tmpMat = zeros(size(tmpMatTime,1)*length(idMetsAll), size(tmpMatTime,2)*length(idMetsAll));
    for m = 1 : length(idMetsAll)
        row = size(tmpMatTime,1)*(m-1) + (1:size(tmpMatTime,1));
        col = size(tmpMatTime,2)*(m-1) + (1:size(tmpMatTime,2));
        tmpMat(row, col) = tmpMat(row, col) + tmpMatTime * tmpMatCoef;
    end
    
    
    row = repmat(length(timeConstr)*((1:length(idMetsAll))-1)', 1, length(tmpTime)) ...
        + repmat((1:length(tmpTime)), length(idMetsAll), 1) ...
        + kk;
    row = reshape(row', numel(row), 1);
        
    colConcConst = nNonPoolMets*(k-1)+idMetsAll;
    colConcRate1 = nNonPoolMets*(nSwitchTimes+1)+nNonPoolMets*(k-1)+idMetsAll;
    colConcRate2 = nNonPoolMets*(nSwitchTimes+1)+nNonPoolMets*(k)+idMetsAll;
    col = [colConcConst, colConcRate1, colConcRate2];
    col = reshape(col', numel(col), 1);
    
    matConcConstRate2Conc(row, col) = matConcConstRate2Conc(row, col) + ...
        1/(fullSwitchTimes(k+1)-fullSwitchTimes(k)) * sparse(tmpMat);
        
    kk = kk + length(tmpTime);
end


end

%% conc -> concMerged
function matConc2ConcMerged = ...
    makeMatConc2ConcMerged(model, expData, optionsMFA, fullSwitchTimes, timeConstr, idMetsAll)

nSwitchTimes = optionsMFA.varSet.nSwitchTimes;

idTimeInterval = zeros(1,length(timeConstr));
switchTimeOrder = 1 : nSwitchTimes+1;
for k = switchTimeOrder
    loc = timeConstr>=fullSwitchTimes(k) & timeConstr<=fullSwitchTimes(k+1);
    idTimeInterval(loc) = k;
end

rowLabelIdMet = zeros(length(timeConstr)*length(idMetsAll),1);
kk = 0;
for k = 1 : nSwitchTimes+1
    loc = idTimeInterval == k;
    tmpTime = timeConstr(loc);
    if isempty(tmpTime)
        continue
    end
    row = repmat(length(timeConstr)*((1:length(idMetsAll))-1)', 1, length(tmpTime)) ...
        + repmat((1:length(tmpTime)), length(idMetsAll), 1) ...
        + kk;
    row = reshape(row', numel(row), 1);
    
    idMetList = repmat(idMetsAll, 1, length(tmpTime));
    idMetList = reshape(idMetList', numel(idMetList), 1);
    
    rowLabelIdMet(row) = idMetList;
    
    kk = kk + length(tmpTime);
end


matConc2ConcMerged = ...
    speye(length(timeConstr)*length(idMetsAll) ,  length(timeConstr)*length(idMetsAll));

idSameEval = optionsMFA.varSet.idSameEval;
unqIdMetSameEval  = unique(idSameEval(idSameEval>=1));
for i = 1 : length(unqIdMetSameEval)
    idMergeMet = find(idSameEval==unqIdMetSameEval(i));
    locMerged = rowLabelIdMet == idMergeMet(1);
    for j = 2 : length(idMergeMet)
        locRemoved = rowLabelIdMet == idMergeMet(j);
        matConc2ConcMerged(locMerged, locRemoved) = ...
            matConc2ConcMerged(locMerged, locRemoved) + speye(nnz(locMerged));
        matConc2ConcMerged(locRemoved,:) = 0; 
    end
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

%% concConst + concRate -> pCoefConc
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

%% param -> concConst
function matParam2ConcConst =  ...
    makeMatParam2ConcConst(model, expData, optionsMFA, fullSwitchTimes)


nNonPoolMets = optionsMFA.varSet.nNonPoolMets;
nSwitchTimes = optionsMFA.varSet.nSwitchTimes;
nParam = optionsMFA.idParamLocal.nParam;

matParam2ConcRate = optionsMFA.convertMat.param2ConcRate;
nConcRates = size(matParam2ConcRate,1); % = nNonPoolMets*(nSwitchTimes+2);
idParamConcRatesMat = reshape((1:nConcRates)', nNonPoolMets, nSwitchTimes+2);

APcell = cell(1,nSwitchTimes+1);
AXcell = cell(1,nSwitchTimes+1);
for k = 1: nSwitchTimes+1
    A_X = sparse(nNonPoolMets, nNonPoolMets);
    if k == 1
        A_P = zeros(nNonPoolMets, nParam);
        row = optionsMFA.idParamLocal.initConcs;
        col = optionsMFA.idParamLocal.initConcs;
        A_P(row, col) = A_P(row, col) + speye(length(row));
    else        
        preA_P = zeros(nNonPoolMets, nConcRates);
        tmpMatTime = [1,  fullSwitchTimes(k),  1/2*fullSwitchTimes(k).^2 ];
        
        tmpMatCoef1 = zeros(3,2);
        tmpMatCoef1(2,1) = fullSwitchTimes(k);
        tmpMatCoef1(2,2) = -fullSwitchTimes(k-1);
        tmpMatCoef1(3,1) = -1;
        tmpMatCoef1(3,2) = 1;
        
        tmpMatCoef2 = zeros(3,2);
        tmpMatCoef2(2,1) = fullSwitchTimes(k+1);
        tmpMatCoef2(2,2) = -fullSwitchTimes(k);
        tmpMatCoef2(3,1) = -1;
        tmpMatCoef2(3,2) = 1;
        
        for m = 1 : nNonPoolMets
            row1 =m;
            col1 =  idParamConcRatesMat(m, k-1:k);
            preA_P(row1, col1) = preA_P(row1, col1) + ...
                1/(fullSwitchTimes(k)-fullSwitchTimes(k-1)) * tmpMatTime * tmpMatCoef1;
            
            row2 = m;
            col2 = idParamConcRatesMat(m, k:k+1);
            preA_P(row2, col2) = preA_P(row2, col2) - ...
                1/(fullSwitchTimes(k+1)-fullSwitchTimes(k)) * tmpMatTime * tmpMatCoef2;
        end
        A_P = preA_P*matParam2ConcRate;

        row = optionsMFA.idParamLocal.initConcs;
        col = optionsMFA.idParamLocal.initConcs;
        A_X(row, col) = A_X(row, col) + speye(length(row));
    end
    APcell{k} = A_P;
    AXcell{k} = A_X;
end

matParam2ConcConst = sparse(nNonPoolMets*(nSwitchTimes+1), nParam);
for k = 1 : nSwitchTimes+1
    tmpAPMat = cat(3, APcell{1:k});
    
    row = nNonPoolMets*(k-1) + (1:nNonPoolMets);
    matParam2ConcConst(row, :) = matParam2ConcConst(row, :) + sum(tmpAPMat, 3);    
end

end

%% 
function [matVar2NotForOpt] = ...
    makeVar2NotForOptim(model, expData, optionsMFA, fullSwitchTimes, varType)

% field2var(optionsMFA.varSet)
idIndFluxes = optionsMFA.varSet.idIndFluxes ;
idNonPoolMets = optionsMFA.varSet.idNonPoolMets;
nRxns = optionsMFA.varSet.nRxns;
switch varType
    case {'indFlux'}
        matNSwitchTimes = optionsMFA.matNSwitchTimesRxns(idIndFluxes,:);
    case {'concRate'}
        matNSwitchTimes = optionsMFA.matNSwitchTimesMets(idNonPoolMets,:);
    case {'nonIndflux'}
        loc = true(nRxns, 1);
        loc(idIndFluxes) = false; 
        matNSwitchTimes = optionsMFA.matNSwitchTimesRxns(loc,:);
    case {'flux'}
        matNSwitchTimes = optionsMFA.matNSwitchTimesRxns;
end

nVar = numel(matNSwitchTimes);
nVarNotForOpt = nnz(~matNSwitchTimes);
matVar2NotForOpt = zeros(nVarNotForOpt, nVar);

rr = 0;
for r = 1 : size(matNSwitchTimes,1)
    if all(matNSwitchTimes(r,:)==true)
        continue
    end
    idFullSwitchTimes = find(~matNSwitchTimes(r,:));
    for k = 1 : length(idFullSwitchTimes)
        idt0 = idFullSwitchTimes(k)-1;
        while true
            t0 = fullSwitchTimes(idt0);
            if matNSwitchTimes(r,idt0)
                break
            end
            idt0 = idt0-1;
        end
        idt2 = idFullSwitchTimes(k)+1;
        while true
            t2 = fullSwitchTimes(idt2);
            if matNSwitchTimes(r,idt2)
                break
            end
            idt2 = idt2+1;
        end
        idt1 = idFullSwitchTimes(k);
        t1 = fullSwitchTimes(idt1);
        
        tmpLoct1 = false(size(matNSwitchTimes));
        tmpLoct1(r,idt1) = true;
        loct1=reshape(tmpLoct1, numel(tmpLoct1),1);
        
        rr = rr +1;
        matVar2NotForOpt(rr, loct1) = 1;
    end
end
matVar2NotForOpt = sparse(matVar2NotForOpt);
end

function [matCalcVar2NotForOpt] = ...
    makeCalcVar2NotForOptim(model, expData, optionsMFA, fullSwitchTimes, varType)

idIndFluxes = optionsMFA.varSet.idIndFluxes ;
idNonPoolMets = optionsMFA.varSet.idNonPoolMets;
nRxns = optionsMFA.varSet.nRxns;
switch varType
    case {'indFlux'}
        matNSwitchTimes = optionsMFA.matNSwitchTimesRxns(idIndFluxes,:);
    case {'concRate'}
        matNSwitchTimes = optionsMFA.matNSwitchTimesMets(idNonPoolMets,:);
    case {'nonIndflux'}
        loc = true(nRxns, 1);
        loc(idIndFluxes) = false; 
        matNSwitchTimes = optionsMFA.matNSwitchTimesRxns(loc,:);
    case {'flux'}
        matNSwitchTimes = optionsMFA.matNSwitchTimesRxns;
end

nVar = numel(matNSwitchTimes);
nVarNotForOpt = nnz(~matNSwitchTimes);
matCalcVar2NotForOpt = zeros(nVarNotForOpt, nVar); 

rr = 0;
for r = 1 : size(matNSwitchTimes,1)
    if all(matNSwitchTimes(r,:)==true)
        continue
    end
    idFullSwitchTimes = find(~matNSwitchTimes(r,:));
    for k = 1 : length(idFullSwitchTimes)
        idt0 = idFullSwitchTimes(k)-1;
        while true
            t0 = fullSwitchTimes(idt0);
            if matNSwitchTimes(r,idt0)
                break
            end
            idt0 = idt0-1;
        end
        idt2 = idFullSwitchTimes(k)+1;
        while true
            t2 = fullSwitchTimes(idt2);
            if matNSwitchTimes(r,idt2)
                break
            end
            idt2 = idt2+1;
        end
        idt1 = idFullSwitchTimes(k);
        t1 = fullSwitchTimes(idt1);
        
        coef = [(t2-t1)/(t2-t0), (t1-t0)/(t2-t0)];
        
        tmpLoct02 = false(size(matNSwitchTimes));
        tmpLoct02(r,[idt0, idt2]) = true;
        loct02=reshape(tmpLoct02, numel(tmpLoct02),1);
        
        rr = rr +1;
        matCalcVar2NotForOpt(rr, loct02) = coef;
    end
end
matCalcVar2NotForOpt = sparse(matCalcVar2NotForOpt);

end

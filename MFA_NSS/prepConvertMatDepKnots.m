%% Make matrix to obtain variales from parameter vector and knot time
function convertMat = prepConvertMatDepKnots171206(model, expData, optionsMFA, fullKnots)
idNonPoolMets = optionsMFA.varSet.idNonPoolMets;
idEvalConcMets = optionsMFA.varSet.idEvalConcMets;
nKnots = optionsMFA.varSet.nKnots;

convertMat = optionsMFA.convertMat;
idParam = optionsMFA.idParamLocal;
nParam = idParam.nParam;
timeConstrQP = linspace(0, expData.time(end), optionsMFA.nStepTimeConstrQP);

idTimeInterval = zeros(1,length(timeConstrQP));
for k = 1 : nKnots+1
    loc = timeConstrQP>=fullKnots(k) & timeConstrQP<=fullKnots(k+1);
    idTimeInterval(loc) = k;    
end

%% param -> concConst
matParam2ConcConst =  ...
    makeMatParam2ConcConst(model, expData, optionsMFA, fullKnots);

%% param -> concRate
matParam2ConcRate=convertMat.param2ConcRate;

%% concConst + concRates -> concs
matConcConstRate2Conc = ...
    makeMatConcConstRate2Conc(model, expData, optionsMFA, fullKnots, timeConstrQP, idNonPoolMets);

%% concConst + concRates -> concsExpTime
matConcConstRate2ConcExpTime = ...
    makeMatConcConstRate2Conc(model, expData, optionsMFA, fullKnots, expData.time, idEvalConcMets);

%% paramLocal -> conc

matParam2Conc = matConcConstRate2Conc * [matParam2ConcConst;matParam2ConcRate];
matParam2ConcExpTime = matConcConstRate2ConcExpTime * [matParam2ConcConst;matParam2ConcRate];


%% param -> pCoefConcs 

matConcConstRate2pCoefConc = ...
    makeMatConcConstRate2pCoefConc(model, expData, optionsMFA, fullKnots, idNonPoolMets);

matParam2pCoefConc = matConcConstRate2pCoefConc * [matParam2ConcConst;matParam2ConcRate];

%% paramLocal -> pCoefFluxes
matFlux2pCoefFlux = ...
    makeMatFlux2pCoefFlux(model, expData, optionsMFA, fullKnots);

matParam2pCoefFlux = matFlux2pCoefFlux * convertMat.param2Flux;

%% Constraint for reaction dependent number of time intervvals
varType = 'concRate';
matConcRate2NotForOptim = ...
    makeVar2NotForOptim(model, expData, optionsMFA, fullKnots, varType);
matCalcConcRate2NotForOptim = ...
    makeCalcVar2NotForOptim(model, expData, optionsMFA, fullKnots, varType);

varType = 'flux';
matFlux2NotForOptim = ...
    makeVar2NotForOptim(model, expData, optionsMFA, fullKnots, varType);
matCalcFlux2NotForOptim = ...
    makeCalcVar2NotForOptim(model, expData, optionsMFA, fullKnots, varType);


matConstrConcRateNotForOptim = ...
    (matCalcConcRate2NotForOptim - matConcRate2NotForOptim) * convertMat.param2ConcRate;
matConstrFluxNotForOptim = ...
    (matCalcFlux2NotForOptim - matFlux2NotForOptim) * convertMat.param2Flux;


%% merge
preConvertMat.param2ConcConst = matParam2ConcConst;
preConvertMat.param2Conc = matParam2Conc;
preConvertMat.param2ConcExpTime = matParam2ConcExpTime;
preConvertMat.concExpTime2ConcMerged = matConcExpTime2ConcMerged;
preConvertMat.param2pCoefConc = matParam2pCoefConc;
preConvertMat.param2pCoefFlux = matParam2pCoefFlux;
preConvertMat.constrConcRateNotForOptim = matConstrConcRateNotForOptim;
preConvertMat.constrFluxNotForOptim = matConstrFluxNotForOptim;

tmpFieldNames = fieldnames(preConvertMat);
for f = 1 : length(tmpFieldNames)
    tmpMat = preConvertMat.(tmpFieldNames{f});
%     disp(tmpFieldNames{f})
%     disp(num2str(nnz(tmpMat)))
%     tmpMat(abs(tmpMat)<=10^-12) = 0;  % 計算誤差を除く⇒この関数では意味がない
%     disp(num2str(nnz(tmpMat)))
    convertMat.(tmpFieldNames{f})= tmpMat;
end

end

%% concConst + concRate -> conc
function matConcConstRate2Conc = ...
    makeMatConcConstRate2Conc(model, expData, optionsMFA, fullKnots, timeConstr, idMetsAll, isEarlierTI)

% field2var(optionsMFA.varSet)
nNonPoolMets = optionsMFA.varSet.nNonPoolMets;
nKnots = optionsMFA.varSet.nKnots;

% if nargin<=6
%     isEarlierTI = false;
% end

idTimeInterval = zeros(1,length(timeConstr));
% if isEarlierTI
%     knotOrder = nKnots+1 : -1 : 1;
% else
    knotOrder = 1 : nKnots+1;
% end
for k = knotOrder
    loc = timeConstr>=fullKnots(k) & timeConstr<=fullKnots(k+1);
    idTimeInterval(loc) = k;
end


% concConstはnNonPoolMets*(nKnots+1), concRatesはnNonPoolMets*(nKnots+2);
matConcConstRate2Conc = sparse(...
    length(timeConstr)*length(idMetsAll),...
    nNonPoolMets*((nKnots+1) + (nKnots+2))...
    ); 
kk = 0;
for k = 1 : nKnots+1 % for loopはmが外、kが中
    tmpMatCoef = zeros(3,3);
    tmpMatCoef(1,1) = fullKnots(k+1)-fullKnots(k);
    tmpMatCoef(2,2) = fullKnots(k+1);
    tmpMatCoef(2,3) = -fullKnots(k);
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
    
    
%     idMet = idMetsAll(m);
%     row = length(timeConstr)*((1:length(idMetsAll))-1)  + kk +  (1:length(tmpTime));
    row = repmat(length(timeConstr)*((1:length(idMetsAll))-1)', 1, length(tmpTime)) ...
        + repmat((1:length(tmpTime)), length(idMetsAll), 1) ...
        + kk;
    row = reshape(row', numel(row), 1);
        
    colConcConst = nNonPoolMets*(k-1)+idMetsAll;
    colConcRate1 = nNonPoolMets*(nKnots+1)+nNonPoolMets*(k-1)+idMetsAll;
    colConcRate2 = nNonPoolMets*(nKnots+1)+nNonPoolMets*(k)+idMetsAll;
    col = [colConcConst, colConcRate1, colConcRate2];
    col = reshape(col', numel(col), 1);
    
    matConcConstRate2Conc(row, col) = matConcConstRate2Conc(row, col) + ...
        1/(fullKnots(k+1)-fullKnots(k)) * sparse(tmpMat);
        
    kk = kk + length(tmpTime);
end

% % matConcConstRate2Concの各行はどの代謝物か？
% rowLabelIdMet = zeros(size(matConcConstRate2Conc,1),1);
% kk = 0;
% for k = 1 : nKnots+1
%     loc = idTimeInterval == k;
%     tmpTime = timeConstr(loc);
%     if isempty(tmpTime)
%         continue
%     end
%     row = repmat(length(timeConstr)*((1:length(idMetsAll))-1)', 1, length(tmpTime)) ...
%         + repmat((1:length(tmpTime)), length(idMetsAll), 1) ...
%         + kk;
%     row = reshape(row', numel(row), 1);
%     
%     idMetList = repmat(idMetsAll, 1, length(tmpTime));
%     idMetList = reshape(idMetList', numel(idMetList), 1);
%     
%     rowLabelIdMet(row) = idMetList;
%     
%     kk = kk + length(tmpTime);
% end

% % 代謝物濃度のmerge(AcCoAなど)
% idSameEval = optionsMFA.varSet.idSameEval;
% unqIdMetSameEval  = unique(idSameEval(idSameEval>=1));
% for i = 1 : length(unqIdMetSameEval)
%     idMergeMet = find(idSameEval==unqIdMetSameEval(i));
%     locMerged = rowLabelIdMet == idMergeMet(1);
%     for j = 2 : length(idMergeMet)
%         locRemoved = rowLabelIdMet == idMergeMet(j);
%         matConcConstRate2Conc(locMerged,:) = ...
%             matConcConstRate2Conc(locMerged,:) + matConcConstRate2Conc(locRemoved,:);
% %         matConcConstRate2Conc(locRemoved,:) = 0; % lbとかに引っかからないようにそのまま残す
%     end
% end

end

%% conc -> concMerged
function matConc2ConcMerged = ...
    makeMatConc2ConcMerged(model, expData, optionsMFA, fullKnots, timeConstr, idMetsAll)

nKnots = optionsMFA.varSet.nKnots;

idTimeInterval = zeros(1,length(timeConstr));
knotOrder = 1 : nKnots+1;
for k = knotOrder
    loc = timeConstr>=fullKnots(k) & timeConstr<=fullKnots(k+1);
    idTimeInterval(loc) = k;
end

% concの各行はどの代謝物か？
rowLabelIdMet = zeros(length(timeConstr)*length(idMetsAll),1);
kk = 0;
for k = 1 : nKnots+1
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

% 代謝物濃度のmerge(AcCoAなど)

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
    makeMatConcConstRate2pCoefConc(model, expData, optionsMFA, fullKnots, idMetsAll)

nNonPoolMets = optionsMFA.varSet.nNonPoolMets;
nKnots = optionsMFA.varSet.nKnots;
idMetsAll = colVec(idMetsAll);

% concConstはnNonPoolMets*(nKnots+1), concRatesはnNonPoolMets*(nKnots+2);
matConcConstRate2pCoefConc = sparse(...
    length(idMetsAll)*3*(nKnots+1),...
    nNonPoolMets*((nKnots+1) + (nKnots+2))...
    ); 
for k = 1 : nKnots+1 % for loopはmが外、kが中
    tmpMatCoef = zeros(3,3);
    tmpMatCoef(1,1) = fullKnots(k+1)-fullKnots(k);
    tmpMatCoef(2,2) = fullKnots(k+1);
    tmpMatCoef(2,3) = -fullKnots(k);
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
%         row = length(idMetsAll)*3*(k-1) + length(idMetsAll)*(0:2) + m;
                
        colConcConst = nNonPoolMets*(k-1)+idMetsAll;
        colConcRate1 = nNonPoolMets*(nKnots+1)+nNonPoolMets*(k-1)+idMetsAll;
        colConcRate2 = nNonPoolMets*(nKnots+1)+nNonPoolMets*(k)+idMetsAll;
        col = [colConcConst, colConcRate1, colConcRate2];
        col = reshape(col', numel(col), 1);
        matConcConstRate2pCoefConc(row, col) = matConcConstRate2pCoefConc(row, col) + ...
            1/(fullKnots(k+1)-fullKnots(k)) * sparse(tmpMat);
        
end

end

%% concConst + concRate -> pCoefConc
function matFlux2pCoefFlux = ...
    makeMatFlux2pCoefFlux(model, expData, optionsMFA, fullKnots)

nRxns = optionsMFA.varSet.nRxns;
nKnots = optionsMFA.varSet.nKnots;

% fluxはnRxns*(nKnots+2);
matFlux2pCoefFlux = sparse(...
    nRxns*2*(nKnots+1),...
    nRxns*(nKnots+2)...
    ); 
for k = 1 : nKnots+1 % for loopはmが外、kが中
    tmpMatCoef = zeros(2,2);
    tmpMatCoef(1,1) = fullKnots(k+1);
    tmpMatCoef(1,2) = -fullKnots(k);
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
%     row = nRxns*2*(k-1) + nRxns*(0:1) + (1:nRxns);
    col = repmat((1:nRxns)', 1, 2) + repmat(nRxns*(k-1:k), nRxns, 1);
    col = reshape(col', numel(col), 1);
    
    matFlux2pCoefFlux(row, col) = matFlux2pCoefFlux(row, col) + ...
        1/(fullKnots(k+1)-fullKnots(k)) * sparse(tmpMat);
end


end

%% param -> concConst
function matParam2ConcConst =  ...
    makeMatParam2ConcConst(model, expData, optionsMFA, fullKnots)
% 時間区間k (k=0,1,2,...)における定数項をXkと書くと、最終的には何らかの行列AとパラメータPを使って
% X_1 = [A_P1, A_X1]*[P, X0]'
% X_2 = [A_P2, A_X2]*[P, X1]'
% X_3 = [A_P3, A_X3]*[P, X2]'
% ...
% と書ける。
% 従って、
% X1=A_P1*P (A_X1は零行列)
% X2 = A_P2*P + A_X2*X1 = (A_P2 + A_X2*A_P1)*P
% X3 = A_P3*P + A_X3*X2 =  (A_P3 + A_X3*A_P2 + A_X3*A_X2*A_P1)*P
% ...
% と書ける。
% 
% 2017/11/8追記
% パラメータベクトル中のconcRatesには、意味のないものが含まれる (NADHとNADの片方など)
% パラメータベクトルを直接使うのではなく、その中のconcRatesのみを取り出して使うべき。

% field2var(optionsMFA.varSet)
nNonPoolMets = optionsMFA.varSet.nNonPoolMets;
nKnots = optionsMFA.varSet.nKnots;
nParam = optionsMFA.idParamLocal.nParam;

matParam2ConcRate = optionsMFA.convertMat.param2ConcRate;
nConcRates = size(matParam2ConcRate,1); % = nNonPoolMets*(nKnots+2);
idParamConcRatesMat = reshape((1:nConcRates)', nNonPoolMets, nKnots+2);

% ある時間区間の定数項とパラメータから次の時点の定数項を算出する行列の作成
APcell = cell(1,nKnots+1);
AXcell = cell(1,nKnots+1);
for k = 1: nKnots+1
    A_X = sparse(nNonPoolMets, nNonPoolMets);
    if k == 1
        A_P = zeros(nNonPoolMets, nParam);
        row = optionsMFA.idParamLocal.initConcs;
        col = optionsMFA.idParamLocal.initConcs;
        A_P(row, col) = A_P(row, col) + speye(length(row));
    else        
        preA_P = zeros(nNonPoolMets, nConcRates);
        tmpMatTime = [1,  fullKnots(k),  1/2*fullKnots(k).^2 ];
        
        tmpMatCoef1 = zeros(3,2);
        tmpMatCoef1(2,1) = fullKnots(k);
        tmpMatCoef1(2,2) = -fullKnots(k-1);
        tmpMatCoef1(3,1) = -1;
        tmpMatCoef1(3,2) = 1;
        
        tmpMatCoef2 = zeros(3,2);
        tmpMatCoef2(2,1) = fullKnots(k+1);
        tmpMatCoef2(2,2) = -fullKnots(k);
        tmpMatCoef2(3,1) = -1;
        tmpMatCoef2(3,2) = 1;
        
        for m = 1 : nNonPoolMets
            row1 =m;
            col1 =  idParamConcRatesMat(m, k-1:k);
            preA_P(row1, col1) = preA_P(row1, col1) + ...
                1/(fullKnots(k)-fullKnots(k-1)) * tmpMatTime * tmpMatCoef1;
            
            row2 = m;
            col2 = idParamConcRatesMat(m, k:k+1);
            preA_P(row2, col2) = preA_P(row2, col2) - ...
                1/(fullKnots(k+1)-fullKnots(k)) * tmpMatTime * tmpMatCoef2;
        end
        A_P = preA_P*matParam2ConcRate;

        row = optionsMFA.idParamLocal.initConcs;
        col = optionsMFA.idParamLocal.initConcs;
        A_X(row, col) = A_X(row, col) + speye(length(row));
    end
    APcell{k} = A_P;
    AXcell{k} = A_X;
end

% パラメータから次の時点の定数項を算出する行列の作成
matParam2ConcConst = sparse(nNonPoolMets*(nKnots+1), nParam);
for k = 1 : nKnots+1
    tmpAPMat = cat(3, APcell{1:k});
    
    row = nNonPoolMets*(k-1) + (1:nNonPoolMets);
    matParam2ConcConst(row, :) = matParam2ConcConst(row, :) + sum(tmpAPMat, 3);    
end

end

%% 変数からForOpt以外のものを取り出す。
function [matVar2NotForOpt] = ...
    makeVar2NotForOptim(model, expData, optionsMFA, fullKnots, varType)

% field2var(optionsMFA.varSet)
idIndFluxes = optionsMFA.varSet.idIndFluxes ;
idNonPoolMets = optionsMFA.varSet.idNonPoolMets;
nRxns = optionsMFA.varSet.nRxns;
switch varType
    case {'indFlux'}
        matNKnots = optionsMFA.matNKnotsRxns(idIndFluxes,:);
    case {'concRate'}
        matNKnots = optionsMFA.matNKnotsMets(idNonPoolMets,:);
    case {'nonIndflux'}
        loc = true(nRxns, 1);
        loc(idIndFluxes) = false; 
        matNKnots = optionsMFA.matNKnotsRxns(loc,:);
    case {'flux'}
        matNKnots = optionsMFA.matNKnotsRxns;
end

nVar = numel(matNKnots);
nVarNotForOpt = nnz(~matNKnots);
% matVar2NotForOpt = sparse(nVarNotForOpt, nVar);
matVar2NotForOpt = zeros(nVarNotForOpt, nVar);

rr = 0;
for r = 1 : size(matNKnots,1)
    if all(matNKnots(r,:)==true)
        continue
    end
    idFullKnots = find(~matNKnots(r,:));
    for k = 1 : length(idFullKnots)
        % t0, t1, t2の同定
        idt0 = idFullKnots(k)-1;
        while true
            t0 = fullKnots(idt0);
            if matNKnots(r,idt0)
                break
            end
            idt0 = idt0-1;
        end
        idt2 = idFullKnots(k)+1;
        while true
            t2 = fullKnots(idt2);
            if matNKnots(r,idt2)
                break
            end
            idt2 = idt2+1;
        end
        idt1 = idFullKnots(k);
        t1 = fullKnots(idt1);
        
%         % 係数の計算
%         coef = [(t2-t1)/(t2-t0), (t1-t0)/(t2-t0)];
        
        % t1の列の場所の同定
        tmpLoct1 = false(size(matNKnots));
        tmpLoct1(r,idt1) = true;
        loct1=reshape(tmpLoct1, numel(tmpLoct1),1);
        
        % 係数の適用
        rr = rr +1;
        matVar2NotForOpt(rr, loct1) = 1;
    end
end
matVar2NotForOpt = sparse(matVar2NotForOpt);
end

%% 変数からForOpt以外のものを計算する。
function [matCalcVar2NotForOpt] = ...
    makeCalcVar2NotForOptim(model, expData, optionsMFA, fullKnots, varType)

% field2var(optionsMFA.varSet)
idIndFluxes = optionsMFA.varSet.idIndFluxes ;
idNonPoolMets = optionsMFA.varSet.idNonPoolMets;
nRxns = optionsMFA.varSet.nRxns;
switch varType
    case {'indFlux'}
        matNKnots = optionsMFA.matNKnotsRxns(idIndFluxes,:);
    case {'concRate'}
        matNKnots = optionsMFA.matNKnotsMets(idNonPoolMets,:);
    case {'nonIndflux'}
        loc = true(nRxns, 1);
        loc(idIndFluxes) = false; 
        matNKnots = optionsMFA.matNKnotsRxns(loc,:);
    case {'flux'}
        matNKnots = optionsMFA.matNKnotsRxns;
end

nVar = numel(matNKnots);
nVarNotForOpt = nnz(~matNKnots);
% matCalcVar2NotForOpt = sparse(nVarNotForOpt, nVar);
matCalcVar2NotForOpt = zeros(nVarNotForOpt, nVar);  % 大きな行列ではないので最後にsparse化

rr = 0;
for r = 1 : size(matNKnots,1)
    if all(matNKnots(r,:)==true)
        continue
    end
    idFullKnots = find(~matNKnots(r,:));
    for k = 1 : length(idFullKnots)
        % t0, t1, t2の同定
        idt0 = idFullKnots(k)-1;
        while true
            t0 = fullKnots(idt0);
            if matNKnots(r,idt0)
                break
            end
            idt0 = idt0-1;
        end
        idt2 = idFullKnots(k)+1;
        while true
            t2 = fullKnots(idt2);
            if matNKnots(r,idt2)
                break
            end
            idt2 = idt2+1;
        end
        idt1 = idFullKnots(k);
        t1 = fullKnots(idt1);
        
        % 係数の計算
        coef = [(t2-t1)/(t2-t0), (t1-t0)/(t2-t0)];
        
        % to, t2の列の場所の同定
        tmpLoct02 = false(size(matNKnots));
        tmpLoct02(r,[idt0, idt2]) = true;
        loct02=reshape(tmpLoct02, numel(tmpLoct02),1);
        
        % 係数の適用
        rr = rr +1;
        matCalcVar2NotForOpt(rr, loct02) = coef;
    end
end
matCalcVar2NotForOpt = sparse(matCalcVar2NotForOpt);

end

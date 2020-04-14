%% Calculation of residuals
function [residuals, fluxEst, factor] = calcResType0(fluxMFA, regLevels, paramUnderOptim, normExpData, paramInfo, optionsKA)

nParam = length(paramUnderOptim);
nData = length(normExpData(1).level);
idRegParam = find(paramInfo.isRegParam);
nReg = length(idRegParam);

paramInfoAll = optionsKA.paramInfo;

%% Parameters
param = zeros(nParam, 1);
for p = 1 : nParam
    tmpFxn = paramInfoAll.fxnGetOriParam{paramInfo.paramRegulatorTypeIds(p)};
    param(p) = tmpFxn(paramUnderOptim(p));
end

loc = strcmp(paramInfoAll.paramTypes(paramInfo.paramTypeIds), 'k1');
k1 = param(loc);
loc = strcmp(paramInfoAll.paramTypes(paramInfo.paramTypeIds), 'Ks');
Ks = param(loc);
loc = strcmp(paramInfoAll.paramTypes(paramInfo.paramTypeIds), 'Kp');
Kp = param(loc);
loc = strcmp(paramInfoAll.paramTypes(paramInfo.paramTypeIds), 'k2');
k2 = param(loc);

if isempty(k1)
    k1 = 0;
end
if isempty(k2)
    k2=0;
end

%% fS
nSubs = length(Ks);
nProd = length(Kp);

if nSubs>=1
    locSubsInParam = paramInfo.paramRegulatorTypeIds==5;
    subsLevel = regLevels(:,locSubsInParam)';
    locSubsInNormExpData = paramInfo.locNormExpData(locSubsInParam);
    subsSCoef = vertcat(normExpData(locSubsInNormExpData).ScoefRegulatorMet);
    S_KsT  = subsLevel ./ repmat(Ks,1,nData);
    if optionsKA.isEstForBack
        S_KsD = repmat(S_KsT(:,1:nData/2),1,2)+1;
    else
        S_KsD = S_KsT+1;
    end
    for i = 1 : nSubs
        S_KsT(i,:) = S_KsT(i,:) .^ subsSCoef(i);
        S_KsD(i,:) = S_KsD(i,:) .^ subsSCoef(i);
    end
else 
    S_KsT= ones(1,nData);
    S_KsD= ones(1,nData);
end
if nProd>=1
    locProdInParam = paramInfo.paramRegulatorTypeIds==6;
    prodLevel = regLevels(:,locProdInParam)';
    locProdInNormExpData = paramInfo.locNormExpData(locProdInParam);
    prodSCoef = vertcat(normExpData(locProdInNormExpData).ScoefRegulatorMet);
    P_KpT = prodLevel ./ repmat(Kp,1,nData);
    if optionsKA.isEstForBack
        P_KpD = repmat(P_KpT(:,nData/2+1:end),1,2) + 1;
    else
        P_KpD = P_KpT + 1;
    end
    for i = 1 : nProd
        P_KpT(i,:) = P_KpT(i,:) .^ prodSCoef(i);
        P_KpD(i,:) = P_KpD(i,:) .^ prodSCoef(i);
    end
else 
    P_KpT= ones(1,nData);  
    P_KpD= ones(1,nData);    
end

Tf =k1*prod(S_KsT,1);
Tb =k2 * prod(P_KpT,1);
D = prod(S_KsD,1) + prod(P_KpD,1) - 1;

fS_f = Tf./D;
fS_b = Tb./D;


%% fR
if nReg == 0
    fR = ones(1,nData);
else
    fR = ones(nReg,nData);
    for r = 1 : nReg
        Kr = param(idRegParam(r));
        tmpRegLevel = regLevels(:,idRegParam(r))';
        if Kr >0 % activation
            fR(r, :) = tmpRegLevel ./ (Kr + tmpRegLevel);
        else % inhibition
            fR(r, :) = -Kr ./ (-Kr + tmpRegLevel);
        end
    end
end

%% flux
if k1>0
    log10AbsFluxEst_f = (log10(abs(fS_f)) + sum(log10(fR),1))';
    fluxEst_f = 10.^log10AbsFluxEst_f;
else
    fluxEst_f = zeros(1,nData)';
end

if k2>0
    log10AbsFluxEst_b = (log10(abs(fS_b)) + sum(log10(fR),1))';
    fluxEst_b = 10.^log10AbsFluxEst_b;
else
    fluxEst_b = zeros(1,nData)';
end

if optionsKA.isEstForBack
    fluxEst_f = fluxEst_f(1:nData/2);
    fluxEst_b = fluxEst_b(nData/2+1:end);
end

factor.fS_f = fS_f;
factor.fS_b = fS_b;
factor.fR = fR;

%% Calculate residuals
if optionsKA.isEstForBack
    fluxEst = [fluxEst_f;fluxEst_b];
else
    fluxEst = [fluxEst_f-fluxEst_b];
end
if all(fluxEst.*fluxMFA>0)
    residuals = log(abs(fluxMFA)) -log(abs(fluxEst));
else
    residuals = zeros(nData,1)+sqrt(optionsKA.optionsOpt.infeasibleScore/nData);
end


end

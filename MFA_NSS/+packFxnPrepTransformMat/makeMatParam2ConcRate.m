%% Make transformation matrix (param -> concRate)
function matParam2ConcRate = ...
    makeMatParam2ConcRate(model, expData, optionsMFA, fullSwitchTimes, inputFlag)

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
        matParam2ConcRate = [];
    case 12
        matParam2ConcRate = sparse(nNonPoolMets*(nSwitchTimes+2), nParamLocal);
        for k = 1 : nSwitchTimes+2
            row = nNonPoolMets*(k-1)+(1:nNonPoolMets);
            col =  idParamLocal.concRates(nNonPoolMets*(k-1)+(1:nNonPoolMets));
            matParam2ConcRate(row, col) = matParam2ConcRate(row, col) + eye(nNonPoolMets);
        end
    case 21
        matParam2ConcRate =  ...
            makeMatParam2ConcRateDepSwitchTimes(model, expData, optionsMFA, fullSwitchTimes);
    case 22
        matParam2ConcRate = optionsMFA.transformMat.param2ConcRate;
end

end

%% param -> concRate (dependent on switch time)
function matParam2ConcRate =  makeMatParam2ConcRateDepSwitchTimes(model, expData, optionsMFA, fullSwitchTimes)

nNonPoolMets = optionsMFA.varSet.nNonPoolMets;
nSwitchTimes = optionsMFA.varSet.nSwitchTimes;
nParam = optionsMFA.idParamLocal.nParam;
idParam = optionsMFA.idParamLocal;
idParamConc = reshape(idParam.concs, nNonPoolMets, nSwitchTimes+2);

Ak = zeros(nSwitchTimes+1,1);
for k = 1 : nSwitchTimes+1
   Ak(k) = 2/diff(fullSwitchTimes(k:k+1)); 
end

matParam2ConcRate = zeros(nNonPoolMets * (nSwitchTimes+2), nParam);
matParam2ConcRate(1:nNonPoolMets, idParam.initConcRates) = eye(nNonPoolMets);
for k = 1 : nSwitchTimes+1
    tmpMat1 = zeros(1, k+1);
    tmpMat1(1) = Ak(1);
    tmpMat1(end) = Ak(k);
    if k+1>=2
        for kk = 2 : k
            tmpMat1(kk) = Ak(kk-1)+Ak(kk);
        end
    end
    
    tmpMat2 = zeros(k+1);
    for kk = 1: k+1
        tmpMat2(kk,kk) = (-1)^(k+1-kk);
    end    
    
    for m = 1 : nNonPoolMets
        row = nNonPoolMets*(k)+m;
        col1 = idParamConc(m,1:k+1);
        matParam2ConcRate(row, col1) = tmpMat1*tmpMat2;
        col2 = idParam.initConcRates(m);
        matParam2ConcRate(row, col2) = (-1)^k;
    end    
    
end

matParam2ConcRate = sparse(matParam2ConcRate);

end
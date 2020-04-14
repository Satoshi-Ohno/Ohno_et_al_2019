%% Make transformation matrix (param -> concConst)
function matParam2ConcConst =  ...
    makeMatParam2ConcConst(model, expData, optionsMFA, fullSwitchTimes, inputFlag);

%% prepare variables
idParamLocal = optionsMFA.idParamLocal;
nParam = idParamLocal.nParam;
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
    case {11, 12}
        matParam2ConcConst = [];
    case {21, 22}
        matParam2ConcRate = optionsMFA.matParam2ConcRate;
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

end
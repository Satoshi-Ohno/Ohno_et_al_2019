%% Make transformation matrix (param -> initConc)
function matParam2InitConc = ...
    makeMatParam2InitConc(model, expData, optionsMFA, fullSwitchTimes, inputFlag)

%% prepare variables
idParamLocal = optionsMFA.idParamLocal;
nParamLocal = idParamLocal.nParam;
nNonPoolMets = optionsMFA.varSet.nNonPoolMets;
idNonPoolMets = optionsMFA.varSet.idNonPoolMets;

%% make matrix
% inputFlag: 
% 11, isempty(fullSwitchTimes) == true & optionsMFA.isUseConcAsParam == true
% 12, isempty(fullSwitchTimes) == true & optionsMFA.isUseConcAsParam == false
% 21, isempty(fullSwitchTimes) == false & optionsMFA.isUseConcAsParam == true
% 22, isempty(fullSwitchTimes) == false & optionsMFA.isUseConcAsParam == false

switch inputFlag
    case {11,12}
        matParam2InitConc  = sparse(nNonPoolMets, nParamLocal);
        row = idNonPoolMets;
        col =  optionsMFA.idParamLocal.initConcs;
        matParam2InitConc(row, col) = matParam2InitConc(row, col) + eye(nNonPoolMets);
    case {21, 22}
        matParam2InitConc = optionsMFA.transformMat.param2InitConc;
end


end

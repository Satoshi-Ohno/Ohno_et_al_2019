%% Add constraints for the same initial metabolite concentrations and fluxes between Ins and Ctrl
function qpModel = modifyQPModelSameInitInsCtrl(model, expData, optionsMFA, qpModel, optimType)
% In the present version, metabolites in the Ins condition must be listed in the first half and metabolites in the Ctrl condition must be listed in the second half.
% The same is true for reactions.
% This will be modified in a future version so that reaction ID and metabolite ID can be used.


nSwitchTimes = optionsMFA.nSwitchTimes;
nNonPoolMets = optionsMFA.varSet.nNonPoolMets;
nRxns = optionsMFA.varSet.nRxns;
nIndFluxes = optionsMFA.varSet.nIndFluxes;

%% Initial metabolite concentrations
idParamInitConcs = optionsMFA.idParamLocal.initConcs;

row = 1:nNonPoolMets/2;
tmpAeqInitConcs = sparse(length(row), optionsMFA.idParamLocal.nParam);
tmpAeqInitConcs(row, idParamInitConcs(1:nNonPoolMets/2)) = speye(length(row));
tmpAeqInitConcs(row, idParamInitConcs(nNonPoolMets/2+1:end)) = -speye(length(row));
tmpbeqInitConcs = zeros(length(row), 1);


%% Initial fluxes
matParam2Flux = optionsMFA.transformMat.param2Flux;
idRowInitFluxes = 1:nRxns;

row = 1:nRxns/2;
preAeqInitFluxes = sparse(length(row), size(matParam2Flux,1));
preAeqInitFluxes(row, idRowInitFluxes(1:nRxns/2)) = speye(length(row));
preAeqInitFluxes(row, idRowInitFluxes(nRxns/2+1:end)) = -speye(length(row));
tmpbeqInitFluxes = zeros(length(row), 1);

tmpAeqInitFluxes = preAeqInitFluxes * matParam2Flux;

%% Add constraints
qpModel.Aeq = [qpModel.Aeq; tmpAeqInitConcs; tmpAeqInitFluxes];
qpModel.beq = [qpModel.beq; tmpbeqInitConcs; tmpbeqInitFluxes];


end

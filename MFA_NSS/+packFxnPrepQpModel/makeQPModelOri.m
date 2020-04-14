%% Make original QP model
function [qpModel] = makeQPModelOri(model, expData, optionsMFA, switchTimes, optimType)
nNonPoolMets = optionsMFA.varSet.nNonPoolMets;
nRxns = optionsMFA.varSet.nRxns;
nSwitchTimes = optionsMFA.varSet.nSwitchTimes;

idParamLocal = optionsMFA.idParamLocal;
nParamLocal = idParamLocal.nParam;
switch optimType
    case {'init', 'metaheuristic'}
        lb = optionsMFA.lbInit;
        ub = optionsMFA.ubInit;
        tmpOptimType = 'init';
    case {'local'}
        lb = optionsMFA.lbLocal;
        ub = optionsMFA.ubLocal;
        tmpOptimType = 'local';
end
fullSwitchTimes = [0, switchTimes', expData.time(end)];

transformMat = optionsMFA.transformMat;

%% Vector of measured metabolite concentrations and the standard deviations
vecExpConcStruct = expData.vecExpConcStruct;

vecExpConcs = vecExpConcStruct.vecExpConcs;
vecSDs = vecExpConcStruct.vecSDs;
vecIsComp = vecExpConcStruct.vecIsComp;
nExpConcs = vecExpConcStruct.nExpConcs;


%% Number of variables
nQPVar = nParamLocal;
idQPVar.paramLocal = 1:nParamLocal;
idQPVar.nVar = nQPVar;

%% Objective function
% RSS=(y-Ax)' * W * (y-Ax)
%       = y'Wy - 2y'WA*x + x'*(A'WA)*x

W = spdiags(1./(vecSDs.^2), 0, nExpConcs, nExpConcs);
W = W *  spdiags(vecIsComp, 0, nExpConcs, nExpConcs);
W(isnan(W)) = 0;

tmpMat = transformMat.param2ConcExpTime;

H = sparse(nQPVar, nQPVar);
H(idQPVar.paramLocal, idQPVar.paramLocal) = ...
    tmpMat' * ...
     W * ...
     tmpMat;
H = H * 2;

q = zeros(1,nQPVar);
q(idQPVar.paramLocal) = -2 *...
    vecExpConcs' * ...
    W * ...
    tmpMat;

%% Constraints for metabolite concentrations
nTmpConstr = size(transformMat.param2ConcQP, 1);
A.concs1 = - transformMat.param2ConcQP;
lbConcs = min(lb.switchTimeConcs,[],2);
tmpLb = diag(lbConcs)*ones(nNonPoolMets, optionsMFA.nStepTimeConstrQP);
b.concs1 = zeros(nTmpConstr, 1) - reshape(tmpLb', nTmpConstr, 1);
ctype.concs1 = repmat('U', nTmpConstr, 1);

A.concs2 =  transformMat.param2ConcQP;
ubConcs = min(ub.switchTimeConcs,[],2);
tmpUb = diag(ubConcs)*ones(nNonPoolMets, optionsMFA.nStepTimeConstrQP);
b.concs2 = zeros(nTmpConstr, 1) + reshape(tmpUb', nTmpConstr, 1);
ctype.concs2 = repmat('U', nTmpConstr, 1);

%% Constraints for change rate in metabolites over time
if optionsMFA.isUseConcAsParam
    nTmpConstr = size(transformMat.param2ConcRate, 1);
    A.concRates1 = -transformMat.param2ConcRate;
    b.concRates1= zeros(nTmpConstr, 1) - reshape(lb.switchTimeConcRates, nNonPoolMets*(nSwitchTimes+2),1);
    A.concRates2 = transformMat.param2ConcRate;
    b.concRates2= zeros(nTmpConstr, 1) + reshape(ub.switchTimeConcRates, nNonPoolMets*(nSwitchTimes+2),1);
    ctype.concRates = repmat('U', nTmpConstr*2, 1);
else
    A.concRates1 = [];
    b.concRates1 = [];
    A.concRates2 = [];
    b.concRates2 = [];
    ctype.concRates = [];
end

%% Constraints for fluxes

nTmpConstr = size(transformMat.param2Flux, 1);
A.flux1 = -transformMat.param2Flux;
b.flux1= zeros(nTmpConstr, 1) - reshape(lb.switchTimeFluxesAll, nRxns*(nSwitchTimes+2),1);
A.flux2 = transformMat.param2Flux;
b.flux2= zeros(nTmpConstr, 1) + reshape(ub.switchTimeFluxesAll, nRxns*(nSwitchTimes+2),1);
ctype.flux = repmat('U', nTmpConstr*2, 1);

%% Bound of variables
lbQP = param2vector(lb, optionsMFA, tmpOptimType);
ubQP = param2vector(ub, optionsMFA, tmpOptimType);

%% Number of constraints
tmpFieldNames= fieldnames(A);
nQPConstr = 0;
for f = 1: length(tmpFieldNames)
    nTmpQPConstr = size(A.(tmpFieldNames{f}),1);
    idQPConstr.(tmpFieldNames{f}) = nQPConstr+(1:nTmpQPConstr);
    nQPConstr = nQPConstr + nTmpQPConstr;    
end
idQPConstr.nConstr = nQPConstr;

idQPEqConstr.nConstr = [];

%% Make QP model
qpModel.H = (H+H')/2;
qpModel.q = q;
qpModel.A =  [A.concs1; A.concs2; A.concRates1; A.concRates2; A.flux1; A.flux2; ];
qpModel.b = [b.concs1; b.concs2; b.concRates1; b.concRates2; b.flux1; b.flux2; ];
qpModel.Aeq = [];
qpModel.beq = [];
qpModel.lb = lbQP;
qpModel.ub = ubQP;
qpModel.idQPVar =idQPVar;
qpModel.idQPConstr = idQPConstr;
qpModel.idQPEqConstr = idQPEqConstr;
qpModel.const = (vecExpConcs)' * W *(vecExpConcs);
qpModel.obj.W = W;
qpModel.obj.matVar2ExpData = transformMat.param2ConcExpTime;
qpModel.obj.expData = vecExpConcs;

end





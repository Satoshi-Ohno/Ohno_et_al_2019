%% Make original QP model
function [qpModel] = makeQPModel171025(model, expData, optionsMFA, switchTimes, optimType)
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

convertMat = optionsMFA.convertMat;

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

tmpMat = convertMat.param2ConcExpTime;

H = sparse(nQPVar, nQPVar);
H(idQPVar.paramLocal, idQPVar.paramLocal) = ...
    tmpMat' * ...
     W * ...
     tmpMat;
H = H * 2;

f = zeros(1,nQPVar);
f(idQPVar.paramLocal) = -2 *...
    vecExpConcs' * ...
    W * ...
    tmpMat;

%% Constraints for metabolite concentrations
nTmpConstr = size(convertMat.param2Conc, 1);
A.concs = - convertMat.param2Conc;
tmpLb = diag(lb.initConcs)*ones(nNonPoolMets, optionsMFA.nStepTimeConstrQP);
b.concs = zeros(nTmpConstr, 1) - reshape(tmpLb', nTmpConstr, 1);
ctype.concs = repmat('U', nTmpConstr, 1);


%% Constraints for fluxes

nTmpConstr = size(convertMat.param2Flux, 1);
A.flux1 = -convertMat.param2Flux;
b.flux1= zeros(nTmpConstr, 1) - reshape(lb.switchTimeFluxesAll, nRxns*(nSwitchTimes+2),1);
A.flux2 = convertMat.param2Flux;
b.flux2= zeros(nTmpConstr, 1) + reshape(ub.switchTimeFluxesAll, nRxns*(nSwitchTimes+2),1);
ctype.flux = repmat('U', nTmpConstr*2, 1);

%% Bound of variables
fxnParam2vector = str2func(optionsMFA.fxnName.param2vector);

lbQP = fxnParam2vector(lb, optionsMFA, tmpOptimType);
ubQP = fxnParam2vector(ub, optionsMFA, tmpOptimType);

%% Number of constraints
nQPConstr = size([A.concs; A.flux1; A.flux2],1);
idQPConstr.concs = 1:size(A.concs,1);
idQPConstr.flux1 = size(A.concs,1) + (1:size(A.flux1,1));
idQPConstr.flux2 = size([A.concs;A.flux1],1) + (1:size(A.flux2,1));
idQPConstr.nConstr = nQPConstr;

idQPEqConstr.nConstr = [];

%% Make QP model
qpModel.H = (H+H')/2;
qpModel.q = f;
qpModel.A =  [A.concs; A.flux1; A.flux2; ];
qpModel.b = [b.concs; b.flux1; b.flux2; ];
qpModel.Aeq = [];
qpModel.beq = [];
qpModel.lb = lbQP;
qpModel.ub = ubQP;
qpModel.idQPVar =idQPVar;
qpModel.idQPConstr = idQPConstr;
qpModel.idQPEqConstr = idQPEqConstr;
qpModel.const = (vecExpConcs)' * W *(vecExpConcs);
qpModel.obj.W = W;
qpModel.obj.matVar2ExpData = convertMat.param2ConcExpTime;
qpModel.obj.expData = vecExpConcs;

end





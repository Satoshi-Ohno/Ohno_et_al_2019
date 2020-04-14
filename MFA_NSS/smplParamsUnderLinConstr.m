%% Sample parameters uniformly under linear constraints
function paramSampled= smplParamsUnderLinConstr(model, expData, optionsMFA, ...
    qpModel, x0)
% initial points are the solution from the QP problem

nSampling =5;
field2var(optionsMFA.varSet)
idParamLocal = optionsMFA.idParamLocal;
nParamLocal = idParamLocal.nParam;

%% Constraints
nParam = length(qpModel.lb);
cprndModel.A = [...
    qpModel.A;...
    -eye(nParam);...
    eye(nParam)];
cprndModel.b = [...
    qpModel.b; ...
    -qpModel.lb; ...
    qpModel.ub];

if ~isempty(qpModel.Aeq)
    cprndModel.A = [...
        cprndModel.A;...
        -qpModel.Aeq;
        qpModel.Aeq];
    cprndModel.b = [...
        cprndModel.b; ...
        -(qpModel.beq-10^-6); ...
        qpModel.beq+10^-6];
end


%% Run cprnd
% options.method = 'hitandrun';
options.method = 'achr';
% options.method = 'gibbs';
options.x0 = x0;
options.discard = size(cprndModel.A,2)*10;
options.discard = 1*10^3;
paramSampled = cprnd(nSampling, cprndModel.A, cprndModel.b, options)';

paramSampled = double(paramSampled);

if isempty(paramSampled)
    paramSampled = [];
    return
end
if any(cprndModel.A*paramSampled(:,nSampling) > cprndModel.b)
    paramSampled = [];
    return
end
paramSampled = paramSampled(:,1);


end

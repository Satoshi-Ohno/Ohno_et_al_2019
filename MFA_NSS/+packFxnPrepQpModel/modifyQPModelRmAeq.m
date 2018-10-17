%% Modify QP model by removing equation constraints
function qpModel = modifyQPModelRmAeq(model, expData, optionsMFA, qpModel, Cx, C0)

if isempty(qpModel.Aeq)
    qpModel.L = Cx;
    qpModel.m = zeros(size(qpModel.A,2),1);
    return
end

%% modify QP model

L = Cx;
m = C0 * qpModel.beq;
q = qpModel.q;
H = qpModel.H;
A = qpModel.A;
b = qpModel.b;

qpModel.H = L' * H * L;
qpModel.H = (qpModel.H+qpModel.H')/2;
qpModel.q = (q*L + m'*H*L);
qpModel.const = qpModel.const + 1/2 * m'*H*m + q*m;

qpModel.A = A*L;
qpModel.b = b - A * m;
isValidA = any(qpModel.A,2);
qpModel.A = qpModel.A(isValidA,:);
qpModel.b = qpModel.b(isValidA);

if ~isempty(qpModel.Aeq)
    qpModel.Aeq = [];
    qpModel.beq = [];
end

qpModel.A = [qpModel.A;-L;L];
qpModel.b = [qpModel.b;-qpModel.lb+m;qpModel.ub-m];
qpModel.lb = [];
qpModel.ub = [];


qpModel.L = L;
qpModel.m = m;


end

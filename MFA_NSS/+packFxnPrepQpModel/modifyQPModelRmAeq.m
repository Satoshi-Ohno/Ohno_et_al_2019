function qpModel = modifyQPModelRmAeq171219(model, expData, optionsMFA, qpModel, Cx, C0)
%% 独立変数の同定

if isempty(qpModel.Aeq)
    qpModel.L = Cx;
    qpModel.m = zeros(size(qpModel.A,2),1);
    return
end

%% qpModelの修正
% x ->Lx + mのとき
% 1/2 * x'Hx   -> 1/2 * x' (L'HL)x + m'HLx + 1/2 * m'Hb;
% q'x  -> q'Lx + q'b
% Ax<b ->   ALx < b-Am
%
% 以上より
% H -> L'HL
% q' -> q'L + m'HL
% const -> const + 1/2 * m'Hm + q'm 
% A -> AL
% b -> b - Am

L = Cx;
m = C0 * qpModel.beq;
% L = tmpMatV;
% m = tmpMatC * tmpBeq;
q = qpModel.q;
H = qpModel.H;
A = qpModel.A;
b = qpModel.b;

% % % テスト
% y = (1:size(L,2))';
% x = L * y+m;
% testy = x(idIndVars);
% % testy = [x(~isVars);x(idVars(idIndVars))];
% testx =  L * y+m;

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


%% 旧パラメータから新パラメータに変換するための行列の作成
% matParam2WOAeq = sparse(length(qpModel.q), length(isVars));
% matParam2WOAeq(1:length(qpModel.q), [find(~isVars),idVars(idIndVars)]) = speye(length(qpModel.q));
% qpModel.matParam2WOAeq = matParam2WOAeq;
% % % テスト
% y = (1:size(L,2))';
% x = L * y+m;
% testy = matParam2WOAeq*x;


% keyboard;


end

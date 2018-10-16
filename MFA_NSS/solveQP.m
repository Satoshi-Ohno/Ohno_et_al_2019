function [qpSol, qpScore, exitflag] = solveQP(qpModel)
options = optimset('Display', 'off', 'Algorithm', 'interior-point-convex');
[qpSol, qpScore, exitflag] = quadprog(...
    qpModel.H, qpModel.q, qpModel.A, qpModel.b ,qpModel.Aeq, qpModel.beq,...
    qpModel.lb, qpModel.ub, [],options);
qpScore = qpScore + qpModel.const;
end
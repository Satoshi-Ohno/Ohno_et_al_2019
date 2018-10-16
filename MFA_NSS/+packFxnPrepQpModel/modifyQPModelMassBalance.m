function qpModel = modifyQPModelMassBalance180116(model, expData, optionsMFA, qpModel, Cx, C0)
%% �������x�̐��񎮂�ǉ�����
% Sv = dC/dt, v = invS * [dC/dt;v^ind]  -> S*invS*[dC/dt;v^ind] - dC/dt = 0
% ���̐��񎮂��Ȃ��ƁANADH��NAD��concRates�̂ǂ������paramLocal�Ɋ܂܂�Ă��܂��B

nConstr = size(optionsMFA.convertMat.constrMassBalance,1);
qpModel.Aeq = [qpModel.Aeq; optionsMFA.convertMat.constrMassBalance];
qpModel.beq = [qpModel.beq; zeros(nConstr,1)];

end

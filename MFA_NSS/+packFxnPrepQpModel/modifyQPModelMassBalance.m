function qpModel = modifyQPModelMassBalance180116(model, expData, optionsMFA, qpModel, Cx, C0)
%% •¨¿ûx‚Ì§–ñ®‚ğ’Ç‰Á‚·‚é
% Sv = dC/dt, v = invS * [dC/dt;v^ind]  -> S*invS*[dC/dt;v^ind] - dC/dt = 0
% ‚±‚Ì§–ñ®‚ª‚È‚¢‚ÆANADH‚ÆNAD‚ÌconcRates‚Ì‚Ç‚¿‚ç‚à‚ªparamLocal‚ÉŠÜ‚Ü‚ê‚Ä‚µ‚Ü‚¤B

nConstr = size(optionsMFA.convertMat.constrMassBalance,1);
qpModel.Aeq = [qpModel.Aeq; optionsMFA.convertMat.constrMassBalance];
qpModel.beq = [qpModel.beq; zeros(nConstr,1)];

end

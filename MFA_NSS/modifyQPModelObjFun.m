%% modify objective function of QP
function qpModel = modifyQPModelObjFun171221(model, expData, optionsMFA, ...
    qpModel, isMinVar)

H = qpModel.H;
nParamQP = length(qpModel.q);
if isMinVar
    qpModel.H=qpModel.H+spdiags(ones(nParamQP,1),0, nParamQP, nParamQP)*10^-6;
else
    if optionsMFA.isUseQPInMH
        qpModel.H=qpModel.H+spdiags(rand(nParamQP,1)*2-1, 0, nParamQP, nParamQP)*10^-6;
    else
        qpModel.H=qpModel.H+spdiags(rand(nParamQP,1)*2-1, 0, nParamQP, nParamQP)*10^0; % 2018.9.13çXêV
    end
end

end
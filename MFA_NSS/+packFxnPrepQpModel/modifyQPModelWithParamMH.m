%% Add constraints for paramMH
function qpModel = modifyQPModelWithParamMH(model, expData, optionsMFA, ...
    qpModel, optimType, paramMH)

convertMat = optionsMFA.convertMat;
switch optimType
    case {'init'}
        idParam = optionsMFA.idParamLocal;
        nParamLocal = idParam.nParam;
        paramMH(convertMat.isLogParamMH) = 10.^paramMH(convertMat.isLogParamMH);
    case {'metaheuristic'} 
        paramMH(convertMat.isLogParamMH) = 10.^paramMH(convertMat.isLogParamMH);
    case {'QP'}
        isLogParamMHInd = convertMat.isLogParamMH(optionsMFA.isIndParams.MH);
        paramMH(isLogParamMHInd) = 10.^paramMH(isLogParamMHInd);
    case {'local'}
        return
end
idParamMH = optionsMFA.idParamMH;


switch optimType
    case {'init'}
        % bound
        tmpAeq1  = sparse(length(idParam.switchTimes),nParamLocal);
        tmpAeq1(1:length(idParam.switchTimes), idParam.switchTimes) = speye(length(idParam.switchTimes));
        tmpBeq1 = paramMH(idParamMH.switchTimes);
        
        qpModel.Aeq = [qpModel.Aeq; tmpAeq1];
        qpModel.beq = [qpModel.beq; tmpBeq1];
        
    case {'metaheuristic', 'QP'}
        tmpMat = convertMat.paramMH2Ind * convertMat.paramLocal2ParamMH;
        tmpAeq = tmpMat;
        tmpBeq = paramMH;
        
        qpModel.Aeq = [qpModel.Aeq; tmpAeq];
        qpModel.beq = [qpModel.beq; tmpBeq];
        
end

end
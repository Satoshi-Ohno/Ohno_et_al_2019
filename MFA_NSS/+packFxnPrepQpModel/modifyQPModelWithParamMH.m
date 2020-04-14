%% Add constraints for paramMH
function qpModel = modifyQPModelWithParamMH(model, expData, optionsMFA, ...
    qpModel, optimType, paramMH)

transformMat = optionsMFA.transformMat;
nNonPoolMets = optionsMFA.varSet.nNonPoolMets;
nSwitchTimes = optionsMFA.varSet.nSwitchTimes;
switch optimType
    case {'init'}
        idParam = optionsMFA.idParamLocal;
        nParamLocal = idParam.nParam;
        paramMH(transformMat.isLogParamMH) = 10.^paramMH(transformMat.isLogParamMH);
    case {'metaheuristic'} 
        paramMH(transformMat.isLogParamMH) = 10.^paramMH(transformMat.isLogParamMH);
    case {'QP'}
        isLogParamMHInd = transformMat.isLogParamMH(optionsMFA.isIndParams.MH);
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
        
        tmpAeq2 = [];
        tmpBeq2 = [];
        if optionsMFA.isUseConcAsParam
            if optionsMFA.isIncludeInitEndConcsInParamMH
                tmpIdParamConcs = reshape(idParam.concs, nNonPoolMets, nSwitchTimes+2);
                tmpIdParamConcs = tmpIdParamConcs(:,[1,nSwitchTimes+2]);
                tmpIdParamMHConcs = reshape(idParamMH.concs, nNonPoolMets, nSwitchTimes+2);
                tmpIdParamMHConcs = tmpIdParamMHConcs(:,[1,nSwitchTimes+2]);                
                tmpAeq2 = sparse(numel(tmpIdParamConcs),nParamLocal);
                tmpAeq2(1:numel(tmpIdParamConcs), tmpIdParamConcs(:)) = speye(numel(tmpIdParamConcs));
                tmpBeq2 = paramMH(tmpIdParamMHConcs(:));
                if optionsMFA.isSameInitInsCtrl
                    tmpAeq2 = tmpAeq2(nNonPoolMets/2+1:end,:);
                    tmpBeq2 = tmpBeq2(nNonPoolMets/2+1:end,:);
                end
            end
        end
        
        qpModel.Aeq = [qpModel.Aeq; tmpAeq1; tmpAeq2];
        qpModel.beq = [qpModel.beq; tmpBeq1; tmpBeq2];
        
        
    case {'metaheuristic', 'QP'}
        tmpMat = transformMat.paramMH2Ind * transformMat.paramLocal2ParamMH;
        tmpAeq = tmpMat;
        tmpBeq = paramMH;
        
        qpModel.Aeq = [qpModel.Aeq; tmpAeq];
        qpModel.beq = [qpModel.beq; tmpBeq];
        
end

end
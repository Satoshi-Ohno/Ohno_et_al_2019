function qpModel = modifyQPModelWithParamMH171206(model, expData, optionsMFA, ...
    qpModel, optimType, paramMH)

convertMat = optionsMFA.convertMat;
% if isForOptim
%     switch optimType
%         case {'init'}
%             idParam = optionsMFA.idParamLocal;
% %             idParam = optionsMFA.idParamLocalForOptim;
% %             paramMH(convertMat.isLogParamMH) = 10.^paramMH(convertMat.isLogParamMH);
%         case {'metaheuristic'}
%             idParam = optionsMFA.idParamMHForOptim;
%     end
%     paramMH(convertMat.isLogParamMHForOptim) = 10.^paramMH(convertMat.isLogParamMHForOptim);
%     idParamMH = optionsMFA.idParamMHForOptim;
% else
    switch optimType
        case {'init'}
            idParam = optionsMFA.idParamLocal;
            nParamLocal = idParam.nParam;
            paramMH(convertMat.isLogParamMH) = 10.^paramMH(convertMat.isLogParamMH);
        case {'metaheuristic'} % MH�ł̓Ɨ��ϐ��̓���Ɏg��
%             idParam = optionsMFA.idParamMH;
            paramMH(convertMat.isLogParamMH) = 10.^paramMH(convertMat.isLogParamMH);
        case {'QP'}
            isLogParamMHInd = convertMat.isLogParamMH(optionsMFA.isIndParams.MH);
            paramMH(isLogParamMHInd) = 10.^paramMH(isLogParamMHInd);
        case {'local'}
            return
    end
    idParamMH = optionsMFA.idParamMH;
% end
% nParam = idParam.nParam;
% paramMH(convertMat.isLogParamMH) = 10.^paramMH(convertMat.isLogParamMH);


switch optimType
    %% �œK�������l����
    % knots��pCoefconcs�Ɋւ���bound���C������B
    % �ړI�֐���noise��������
    case {'init'}
        % bound
        tmpAeq1  = sparse(length(idParam.knots),nParamLocal);
        tmpAeq1(1:length(idParam.knots), idParam.knots) = speye(length(idParam.knots));
        tmpBeq1 = paramMH(idParamMH.knots);
        
        tmpAeq2  = zeros(length(idParam.coefCorrCompt),nParamLocal);
        tmpAeq2(1:length(idParam.coefCorrCompt), idParam.coefCorrCompt) = speye(length(idParam.coefCorrCompt));
        tmpBeq2 = paramMH(idParamMH.coefCorrCompt);
        
        qpModel.Aeq = [qpModel.Aeq; tmpAeq1; tmpAeq2];
        qpModel.beq = [qpModel.beq; tmpBeq1; tmpBeq2];
        
    %% Metaheuristic�œK��
    % paramMH�͒萔�Ƃ݂Ȃ��AparamLocal���̕ʃp�����[�^��QP�ŋ��߂����B
    % �ړI�֐�, �������, bound���C������
    case {'metaheuristic', 'QP'}
%         if isForOptim
%             tmpMat = ...
%                 convertMat.paramMH2ForOptim* ...
%                 convertMat.paramLocal2ParamMH ;
%         else
            tmpMat = convertMat.paramMH2Ind * convertMat.paramLocal2ParamMH;
%         end
        tmpAeq = tmpMat;
        tmpBeq = paramMH;
        
        qpModel.Aeq = [qpModel.Aeq; tmpAeq];
        qpModel.beq = [qpModel.beq; tmpBeq];
        
end

end
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
        case {'metaheuristic'} % MHでの独立変数の同定に使う
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
    %% 最適化初期値推定
    % knotsとpCoefconcsに関するboundを修正する。
    % 目的関数にnoiseを加える
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
        
    %% Metaheuristic最適化
    % paramMHは定数とみなし、paramLocal中の別パラメータをQPで求めたい。
    % 目的関数, 制約条件, boundを修正する
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
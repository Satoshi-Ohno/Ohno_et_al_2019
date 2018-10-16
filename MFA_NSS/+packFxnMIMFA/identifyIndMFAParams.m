function [isIndMFAParams, isRedundantConstrAeq] = ...
    identifyIndMFAParams171225(model, expData, optionsMFA, optimType)

field2var(optionsMFA.varSet)
idParamLocal = optionsMFA.idParamLocal;
nParamLocal = idParamLocal.nParam;
lbInit = optionsMFA.lbInit;
ubInit = optionsMFA.ubInit;
optimTypeInit = 'init';
idParamMH = optionsMFA.idParamMH;
nParamMH = idParamMH.nParam;

%% paramVectorの準備
tmpParamMH = rand(nParamMH, 1)*3;

%% knotsとcoefCorrComptの初期値
% knot time
log10Range = log10(ubInit.knots)-log10(lbInit.knots);
while true
    knots = 10.^(log10Range.*rand(nKnots, 1)+log10(lbInit.knots));
    knots = sort(knots);
    fullKnots = [0, knots', expData.time(end)];
    if all(diff(fullKnots) >= optionsMFA.minKnotTimeDiff)
        break
    end
end

% compartment間をつなぐ係数
if ~isempty(idComptCorrParam)
    log10Range = log10(ubInit.coefCorrCompt)-log10(lbInit.coefCorrCompt);
    coefCorrCompt = 10.^(log10Range.*rand(nComptCorrParam, 1)+log10(lbInit.coefCorrCompt));
else
    if optionsMFA.isAccountMediaDrop
        coefCorrCompt = []; % media dropを考慮する場合は時間によって補正値が異なる。
    else
        mediaInfo = optionsMFA.mediaInfo;
        coefCorrCompt = 1/mediaInfo.mediaVol*mediaInfo.cellAmount;
    end
end

tmpParamMH(idParamMH.knots) = log10(knots);
if ~isempty(idParamMH.coefCorrCompt)
    tmpParamMH(idParamMH.coefCorrCompt) = log10(coefCorrCompt);
end

%% paramMHを独立パラメータに変換
switch optimType
    case {'QP'}
        tmpParamMH = optionsMFA.convertMat.paramMH2Ind * tmpParamMH;
end

%% qpModelの作成
switch optimType
    case {'QP'}
        isInputIndVars = true;
    otherwise
        isInputIndVars = false;
end
isInputFullIndVarsConstr= false;
fxnPrepQpModel = str2func(optionsMFA.fxnName.prepQpModel);
[qpModel, ~, qpModelConstr] = fxnPrepQpModel(model, expData, optionsMFA,  ...
    optimType, tmpParamMH,isInputIndVars, isInputFullIndVarsConstr);

%% 目的関数の修正
% initではランダムなノイズを加える
% metaheuristicではパラメータの二乗和の最小化
isMinVar = true;
fxnModifyQPModelObjFun = str2func(optionsMFA.fxnName.modifyQPModelObjFun);
qpModel = fxnModifyQPModelObjFun(model, expData, optionsMFA, qpModel, isMinVar);
%     qpModelConstr = fxnModifyQPModelObjFun(model, expData, optionsMFA, qpModelConstr, isMinVar);

%% 代謝物濃度に対してQPでフィッティング
[qpSol, qpScore,exitFlag] = solveQP(qpModel);
% keyboard;
%     [qpSolConstr, qpScoreConstr, exitflag] = solveQP(qpModelConstr);
%     if ~isempty(qpSol)
%     end

isIndMFAParams = false(nParamLocal,1);
isIndMFAParams(qpModel.idIndVars) = true;
switch optimType
    case {'init', 'local', 'QP'}
    case {'metaheuristic'}
        isIndMFAParams = optionsMFA.convertMat.paramLocal2ParamMH * isIndMFAParams==true;
end
isRedundantConstrAeq = false(size(qpModelConstr.Aeq,1), 1);
isRedundantConstrAeq(qpModel.idRedundantConstr) = true;

%% 独立制約式のチェック
% QPでの独立の制約式は、rxnSpecificNKnotsに関するものを優先する。

end

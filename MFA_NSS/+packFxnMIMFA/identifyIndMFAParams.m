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

%% paramVector�̏���
tmpParamMH = rand(nParamMH, 1)*3;

%% knots��coefCorrCompt�̏����l
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

% compartment�Ԃ��Ȃ��W��
if ~isempty(idComptCorrParam)
    log10Range = log10(ubInit.coefCorrCompt)-log10(lbInit.coefCorrCompt);
    coefCorrCompt = 10.^(log10Range.*rand(nComptCorrParam, 1)+log10(lbInit.coefCorrCompt));
else
    if optionsMFA.isAccountMediaDrop
        coefCorrCompt = []; % media drop���l������ꍇ�͎��Ԃɂ���ĕ␳�l���قȂ�B
    else
        mediaInfo = optionsMFA.mediaInfo;
        coefCorrCompt = 1/mediaInfo.mediaVol*mediaInfo.cellAmount;
    end
end

tmpParamMH(idParamMH.knots) = log10(knots);
if ~isempty(idParamMH.coefCorrCompt)
    tmpParamMH(idParamMH.coefCorrCompt) = log10(coefCorrCompt);
end

%% paramMH��Ɨ��p�����[�^�ɕϊ�
switch optimType
    case {'QP'}
        tmpParamMH = optionsMFA.convertMat.paramMH2Ind * tmpParamMH;
end

%% qpModel�̍쐬
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

%% �ړI�֐��̏C��
% init�ł̓����_���ȃm�C�Y��������
% metaheuristic�ł̓p�����[�^�̓��a�̍ŏ���
isMinVar = true;
fxnModifyQPModelObjFun = str2func(optionsMFA.fxnName.modifyQPModelObjFun);
qpModel = fxnModifyQPModelObjFun(model, expData, optionsMFA, qpModel, isMinVar);
%     qpModelConstr = fxnModifyQPModelObjFun(model, expData, optionsMFA, qpModelConstr, isMinVar);

%% ��ӕ��Z�x�ɑ΂���QP�Ńt�B�b�e�B���O
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

%% �Ɨ����񎮂̃`�F�b�N
% QP�ł̓Ɨ��̐��񎮂́ArxnSpecificNKnots�Ɋւ�����̂�D�悷��B

end

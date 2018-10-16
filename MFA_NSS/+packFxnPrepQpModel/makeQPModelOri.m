function [qpModel] = makeQPModel171025(model, expData, optionsMFA, knots, coefCorrCompt, optimType)
%%  ����
% 2017.10.23 concConst�̕ϐ��͗p�ӂ��Ȃ��B
% 2017.10.24 paramMH�Ō��܂��Ă���Ƃ���͌Œ肷��
% field2var(optionsMFA.varSet)
nNonPoolMets = optionsMFA.varSet.nNonPoolMets;
nRxns = optionsMFA.varSet.nRxns;
nKnots = optionsMFA.varSet.nKnots;

idParamLocal = optionsMFA.idParamLocal;
nParamLocal = idParamLocal.nParam;
switch optimType
    case {'init', 'metaheuristic'}
        lb = optionsMFA.lbInit;
        ub = optionsMFA.ubInit;
        tmpOptimType = 'init';
    case {'local'}
        lb = optionsMFA.lbLocal;
        ub = optionsMFA.ubLocal;
        tmpOptimType = 'local';
end
fullKnots = [0, knots', expData.time(end)];

convertMat = optionsMFA.convertMat;

%% QP�ő�ӕ��Z�x�ɂ��Ẵt�B�b�e�B���O
% �ړI�֐��F�����ӕ���(���U�ŏd�ݕt������)�c�����a
% �ϐ��FparamLocal�ƃV�~�����[�V�����̑�ӕ��Z�x�̃x�N�g��
% �������: concs, knotConcRates, flux��bound��

%% �����f�[�^�ƕ��U�̃x�N�g��
vecExpConcStruct = expData.vecExpConcStruct;
vecExpConcStruct = modifyVecExpConcStructWithCoefCorrCompt(...
    model, expData, optionsMFA, vecExpConcStruct, coefCorrCompt);

vecExpConcs = vecExpConcStruct.vecExpConcs;
vecSDs = vecExpConcStruct.vecSDs;
vecIsComp = vecExpConcStruct.vecIsComp;
nExpConcs = vecExpConcStruct.nExpConcs;


%% �ϐ��̐�
% local�̃p�����[�^
% nKnots�̐������A��ӕ��Z�x�̒萔�� -> ��ӕ��Z�x�̎Z�o
% �e���ԋ�Ԃɂ��āAknot���Ԃł̑�ӕ��Z�x -> ���� ->���񎮂̂ݎc���B�ϐ��͎c���Ȃ��B
% �������Ԃ̑�ӕ��Z�x -> �c���̌v�Z

nQPVar = nParamLocal;
idQPVar.paramLocal = 1:nParamLocal;
idQPVar.nVar = nQPVar;

%% �ړI�֐� (�c���v�Z�p�̑�ӕ��Z�x�̕ϐ����g��Ȃ��ꍇ)
% y: �����f�[�^, x: �p�����[�^, A: �p�����[�^-->��ӕ��Z�x�̕ϊ��s��, W�F�c���̏d�ݕt��
% RSS=(y-Ax)' * W * (y-Ax)
%       = y'Wy - 2y'WA*x + x'*(A'WA)*x
% ��1���͒萔, ��2����x��1���̍��A��3����x��2���̍�
% 

W = spdiags(1./(vecSDs.^2), 0, nExpConcs, nExpConcs);
W = W *  spdiags(vecIsComp, 0, nExpConcs, nExpConcs);
W(isnan(W)) = 0;

tmpMat = convertMat.concExpTime2ConcMerged * convertMat.param2ConcExpTime;

H = sparse(nQPVar, nQPVar);
H(idQPVar.paramLocal, idQPVar.paramLocal) = ...
    tmpMat' * ...
     W * ...
     tmpMat;
H = H * 2;

f = zeros(1,nQPVar);
f(idQPVar.paramLocal) = -2 *...
    vecExpConcs' * ...
    W * ...
    tmpMat;
% f2(isnan(f2))=0;


%% concs�̐���
nTmpConstr = size(convertMat.param2Conc, 1);
A.concs = - convertMat.param2Conc;
tmpLb = diag(lb.initConcs)*ones(nNonPoolMets, optionsMFA.nStepTimeConstrQP);
b.concs = zeros(nTmpConstr, 1) - reshape(tmpLb', nTmpConstr, 1);
ctype.concs = repmat('U', nTmpConstr, 1);


%% Flux�̐���

nTmpConstr = size(convertMat.param2Flux, 1);
A.flux1 = -convertMat.param2Flux;
b.flux1= zeros(nTmpConstr, 1) - reshape(lb.knotFluxesAll, nRxns*(nKnots+2),1);
% b.flux1= zeros(nTmpConstr, 1) - reshape(optionsMFA.lbLocal.knotFluxesAll, nRxns*(nKnots+2),1);
A.flux2 = convertMat.param2Flux;
b.flux2= zeros(nTmpConstr, 1) + reshape(ub.knotFluxesAll, nRxns*(nKnots+2),1);
% b.flux2= zeros(nTmpConstr, 1) + reshape(optionsMFA.ubLocal.knotFluxesAll, nRxns*(nKnots+2),1);
ctype.flux = repmat('U', nTmpConstr*2, 1);

%% �������L��nKnots�ɔ����������
% ����͕ʂ̊֐��ōs���B
Aeq.rxnSpecificNKnots = [];
beq.rxnSpecificNKnots = [];

%% �ϐ���bound
fxnParam2vector = str2func(optionsMFA.fxnName.param2vector);

%param2vector�ł�optimType��metaheuristic�ł͂��߁Blocal�̒����̃x�N�g���ɂȂ�Ȃ��B
lbQP = fxnParam2vector(lb, optionsMFA, tmpOptimType);
ubQP = fxnParam2vector(ub, optionsMFA, tmpOptimType);

%% ���񎮂̐�
nQPConstr = size([A.concs; A.flux1; A.flux2],1);
idQPConstr.concs = 1:size(A.concs,1);
idQPConstr.flux1 = size(A.concs,1) + (1:size(A.flux1,1));
idQPConstr.flux2 = size([A.concs;A.flux1],1) + (1:size(A.flux2,1));
% nQPConstr = size([A.concs; A.flux1; A.flux2; A.rxnSpecificNKnots1; A.rxnSpecificNKnots2],1);
% idQPConstr.concs = 1:size(A.concs,1);
% idQPConstr.flux1 = size(A.concs,1) + (1:size(A.flux1,1));
% idQPConstr.flux2 = size([A.concs;A.flux1],1) + (1:size(A.flux2,1));
% idQPConstr.rxnSpecificNKnots1 = size([A.concs;A.flux1;A.flux2],1) + (1:size(A.rxnSpecificNKnots1,1));
% idQPConstr.rxnSpecificNKnots2 = size([A.concs;A.flux1;A.flux2;A.rxnSpecificNKnots1],1) + (1:size(A.rxnSpecificNKnots2,1));
idQPConstr.nConstr = nQPConstr;

% nQPEqConstr = size(Aeq.rxnSpecificNKnots,1);
% idQPEqConstr.rxnSpecificNKnots= 1:size(Aeq.rxnSpecificNKnots,1);
idQPEqConstr.nConstr = [];

%% QP model�̍쐬
% quadprog
qpModel.H = (H+H')/2;
qpModel.q = f;
% qpModel.A =  [A.flux1; A.flux2];
% qpModel.b = [b.flux1; b.flux2];
qpModel.A =  [A.concs; A.flux1; A.flux2; ];
qpModel.b = [b.concs; b.flux1; b.flux2; ];
% qpModel.A =  [A.concs; A.flux1; A.flux2; A.rxnSpecificNKnots1; A.rxnSpecificNKnots2];
% qpModel.b = [b.concs; b.flux1; b.flux2;  b.rxnSpecificNKnots1; b.rxnSpecificNKnots2];
% qpModel.Aeq = [];
% qpModel.beq = [];
qpModel.Aeq = [Aeq.rxnSpecificNKnots];
qpModel.beq = [beq.rxnSpecificNKnots];
qpModel.lb = lbQP;
qpModel.ub = ubQP;
qpModel.idQPVar =idQPVar;
qpModel.idQPConstr = idQPConstr;
qpModel.idQPEqConstr = idQPEqConstr;
qpModel.const = (vecExpConcs)' * W *(vecExpConcs);
qpModel.obj.W = W;
qpModel.obj.matVar2ExpData = convertMat.param2ConcExpTime;
qpModel.obj.expData = vecExpConcs;


% options = optimset('Display', 'off');
% [qpSol, qpScore, exitflag] = quadprog(...
%     qpModel.H, qpModel.q, qpModel.A, qpModel.b ,qpModel.Aeq, qpModel.beq,...
%     qpModel.lb, qpModel.ub, [],options);
% qpScore = qpScore + qpModel.const;
% keyboard;
% % qpScore = qpScore+ sum(vecExpConcs(vecSDs~=0).^2./(vecSDs(vecSDs~=0).^2));
% qpScore = qpScore  + (vecExpConcs)' * W *(vecExpConcs);

% testScore = (vecExpConcs)' * W *(vecExpConcs) +...
%     f'*qpSol + 0.5 * qpSol' * H * qpSol;
% 
% ansScore = (vecExpConcs-convertMat.param2ConcExpTime*qpSol)' ...
%     * W *...
%     (vecExpConcs-convertMat.param2ConcExpTime*qpSol);
% keyboard;

end

%% ��ӕ��Z�x�̂���SD���R���p�[�g�����g�␳�l�ɏ]���ĕ␳
function vecExpConcStruct = modifyVecExpConcStructWithCoefCorrCompt(...
    model, expData, optionsMFA, vecExpConcStruct, coefCorrCompt)

% field2var(optionsMFA.varSet)
% field2var(vecExpConcStruct);
isEmptyCoefCorrCompt = isempty(coefCorrCompt);

for m = 1 : vecExpConcStruct.nExpConcs
    idMet = vecExpConcStruct.vecIdMet(m);
    tmpExpConc = vecExpConcStruct.vecExpConcs(m);
    tmpSD = vecExpConcStruct.vecSDs(m);
    
    % compartment�ɂ��␳
    [tmpExpConc, tmpSD] = corrComptMedia2Cell(model, expData, optionsMFA, ...
        tmpExpConc, tmpSD, coefCorrCompt, idMet, isEmptyCoefCorrCompt);

    % �܂Ƃ�
    vecExpConcStruct.vecExpConcs(m) = tmpExpConc;
    vecExpConcStruct.vecSDs(m) = tmpSD;
end

end


%% ��ӕ��Z�x�̂���SD���R���p�[�g�����g�␳�l�ɏ]���ĕ␳
function     [tmpExpConcs, tmpSDs] = corrComptMedia2Cell(model, expData, optionsMFA, ...
    tmpExpConcs, tmpSDs, paramCoefCorrCompt, idMet, isEmptyCoefCorrCompt)

% field2var(optionsMFA.varSet)
nExpTime = optionsMFA.varSet.nExpTime;

idCompt = model.metInfo.compt(idMet);
if isEmptyCoefCorrCompt
    if idCompt>=2
        field2var(optionsMFA.mediaInfo)
        isComp = expData.concs(idMet,:)>0;  % 0����thresholdMDV�ɏC��
        
        % �|�n�Z�x�̌o���ω�
        tmpIdExpTime = (0:nExpTime-1);
%         tmpIdExpTime = tmpIdExpTime(isComp);
        mediaVolSeries = mediaVol-tmpIdExpTime*mediaDropPerSampling;
        
        % ��ӕ��Z�x�̕␳
        amount = tmpExpConcs .* mediaVolSeries'; % nmol
        smplAmount = [0;tmpExpConcs(1:end-1)]*mediaDropPerSampling; % �T���v�����O�ŏ����ꂽnmol
        smplAmount(isnan(smplAmount)) = 0;
        tmpExpConcs = (amount+cumsum(smplAmount)) / cellAmount;
        
        % ��ӕ��Z�x��SD�̕␳
        amountVar = tmpSDs .* mediaVolSeries'; % nmol
        smplAmountSD = [0;tmpSDs(1:end-1)]*(mediaDropPerSampling); % �T���v�����O�ŏ����ꂽnmol
        smplAmountSD(isnan(smplAmountSD)) = 0;
        tmpSDs= (amountVar+cumsum(smplAmountSD)) / cellAmount;
        
    else
        coefCorrCompt = [1;paramCoefCorrCompt];
        tmpExpConcs = tmpExpConcs / coefCorrCompt(idCompt);
        tmpSDs = tmpSDs /  coefCorrCompt(idCompt);
    end
else
    coefCorrCompt = [1;paramCoefCorrCompt];
    tmpExpConcs = tmpExpConcs / coefCorrCompt(idCompt);
    tmpSDs = tmpSDs /  coefCorrCompt(idCompt);
end

end





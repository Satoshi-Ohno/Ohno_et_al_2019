function [qpModel] = makeQPModel171025(model, expData, optionsMFA, knots, coefCorrCompt, optimType)
%%  準備
% 2017.10.23 concConstの変数は用意しない。
% 2017.10.24 paramMHで決まっているところは固定する
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

%% QPで代謝物濃度についてのフィッティング
% 目的関数：測定代謝物の(分散で重み付けした)残差二乗和
% 変数：paramLocalとシミュレーションの代謝物濃度のベクトル
% 制約条件: concs, knotConcRates, fluxがbound内

%% 実験データと分散のベクトル
vecExpConcStruct = expData.vecExpConcStruct;
vecExpConcStruct = modifyVecExpConcStructWithCoefCorrCompt(...
    model, expData, optionsMFA, vecExpConcStruct, coefCorrCompt);

vecExpConcs = vecExpConcStruct.vecExpConcs;
vecSDs = vecExpConcStruct.vecSDs;
vecIsComp = vecExpConcStruct.vecIsComp;
nExpConcs = vecExpConcStruct.nExpConcs;


%% 変数の数
% localのパラメータ
% nKnotsの数だけ、代謝物濃度の定数項 -> 代謝物濃度の算出
% 各時間区間について、knot時間での代謝物濃度 -> 制約式 ->制約式のみ残す。変数は残さない。
% 実験時間の代謝物濃度 -> 残差の計算

nQPVar = nParamLocal;
idQPVar.paramLocal = 1:nParamLocal;
idQPVar.nVar = nQPVar;

%% 目的関数 (残差計算用の代謝物濃度の変数を使わない場合)
% y: 実験データ, x: パラメータ, A: パラメータ-->代謝物濃度の変換行列, W：残差の重み付け
% RSS=(y-Ax)' * W * (y-Ax)
%       = y'Wy - 2y'WA*x + x'*(A'WA)*x
% 第1項は定数, 第2項はxの1次の項、第3項はxの2次の項
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


%% concsの制約式
nTmpConstr = size(convertMat.param2Conc, 1);
A.concs = - convertMat.param2Conc;
tmpLb = diag(lb.initConcs)*ones(nNonPoolMets, optionsMFA.nStepTimeConstrQP);
b.concs = zeros(nTmpConstr, 1) - reshape(tmpLb', nTmpConstr, 1);
ctype.concs = repmat('U', nTmpConstr, 1);


%% Fluxの制約式

nTmpConstr = size(convertMat.param2Flux, 1);
A.flux1 = -convertMat.param2Flux;
b.flux1= zeros(nTmpConstr, 1) - reshape(lb.knotFluxesAll, nRxns*(nKnots+2),1);
% b.flux1= zeros(nTmpConstr, 1) - reshape(optionsMFA.lbLocal.knotFluxesAll, nRxns*(nKnots+2),1);
A.flux2 = convertMat.param2Flux;
b.flux2= zeros(nTmpConstr, 1) + reshape(ub.knotFluxesAll, nRxns*(nKnots+2),1);
% b.flux2= zeros(nTmpConstr, 1) + reshape(optionsMFA.ubLocal.knotFluxesAll, nRxns*(nKnots+2),1);
ctype.flux = repmat('U', nTmpConstr*2, 1);

%% 反応特有のnKnotsに伴う制約条件
% これは別の関数で行う。
Aeq.rxnSpecificNKnots = [];
beq.rxnSpecificNKnots = [];

%% 変数のbound
fxnParam2vector = str2func(optionsMFA.fxnName.param2vector);

%param2vectorでのoptimTypeがmetaheuristicではだめ。localの長さのベクトルにならない。
lbQP = fxnParam2vector(lb, optionsMFA, tmpOptimType);
ubQP = fxnParam2vector(ub, optionsMFA, tmpOptimType);

%% 制約式の数
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

%% QP modelの作成
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

%% 代謝物濃度のそのSDをコンパートメント補正値に従って補正
function vecExpConcStruct = modifyVecExpConcStructWithCoefCorrCompt(...
    model, expData, optionsMFA, vecExpConcStruct, coefCorrCompt)

% field2var(optionsMFA.varSet)
% field2var(vecExpConcStruct);
isEmptyCoefCorrCompt = isempty(coefCorrCompt);

for m = 1 : vecExpConcStruct.nExpConcs
    idMet = vecExpConcStruct.vecIdMet(m);
    tmpExpConc = vecExpConcStruct.vecExpConcs(m);
    tmpSD = vecExpConcStruct.vecSDs(m);
    
    % compartmentによる補正
    [tmpExpConc, tmpSD] = corrComptMedia2Cell(model, expData, optionsMFA, ...
        tmpExpConc, tmpSD, coefCorrCompt, idMet, isEmptyCoefCorrCompt);

    % まとめ
    vecExpConcStruct.vecExpConcs(m) = tmpExpConc;
    vecExpConcStruct.vecSDs(m) = tmpSD;
end

end


%% 代謝物濃度のそのSDをコンパートメント補正値に従って補正
function     [tmpExpConcs, tmpSDs] = corrComptMedia2Cell(model, expData, optionsMFA, ...
    tmpExpConcs, tmpSDs, paramCoefCorrCompt, idMet, isEmptyCoefCorrCompt)

% field2var(optionsMFA.varSet)
nExpTime = optionsMFA.varSet.nExpTime;

idCompt = model.metInfo.compt(idMet);
if isEmptyCoefCorrCompt
    if idCompt>=2
        field2var(optionsMFA.mediaInfo)
        isComp = expData.concs(idMet,:)>0;  % 0からthresholdMDVに修正
        
        % 培地濃度の経時変化
        tmpIdExpTime = (0:nExpTime-1);
%         tmpIdExpTime = tmpIdExpTime(isComp);
        mediaVolSeries = mediaVol-tmpIdExpTime*mediaDropPerSampling;
        
        % 代謝物濃度の補正
        amount = tmpExpConcs .* mediaVolSeries'; % nmol
        smplAmount = [0;tmpExpConcs(1:end-1)]*mediaDropPerSampling; % サンプリングで除かれたnmol
        smplAmount(isnan(smplAmount)) = 0;
        tmpExpConcs = (amount+cumsum(smplAmount)) / cellAmount;
        
        % 代謝物濃度のSDの補正
        amountVar = tmpSDs .* mediaVolSeries'; % nmol
        smplAmountSD = [0;tmpSDs(1:end-1)]*(mediaDropPerSampling); % サンプリングで除かれたnmol
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





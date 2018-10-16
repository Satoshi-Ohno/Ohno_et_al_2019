function inputData = loadMFAInput(arrayJobId)
%MFAの入力

% 1条件あたりの最適化の試行回数
optionsMFA.nOptimIter = 10; % arrayJobIdでも指定

% internal knots数の種類数
typeNKnots = 0:2;
% typeNKnots = 2;
% typeNKnots = -1;  %反応ごとにnKnotsを指定する
nTypeNKnots = length(typeNKnots);

%実験データの種類
expTypeAll = {'Ins', 'Ctrl'};
nExpTypeAll = length(expTypeAll);

%% arrayJobIDからinternal knotsの数, expTypeの指定
arrayJobIdsAll = 1:optionsMFA.nOptimIter*nTypeNKnots*nExpTypeAll;
matArrayJobIdsAll = reshape(arrayJobIdsAll, optionsMFA.nOptimIter, nTypeNKnots, nExpTypeAll);

for i = 1:nExpTypeAll
    if any(any(matArrayJobIdsAll(:,:,i)==arrayJobId))        
        idExpType = i;
        [idOptimIter, idTypeNKnots] = find(matArrayJobIdsAll(:,:,i)==arrayJobId);
        break
    end
end

% internal knotsの数の指定
optionsMFA.nKnots = []; 
% arrayJobIdでnKnotsを修正
if isempty(optionsMFA.nKnots)
    optionsMFA.nKnots = typeNKnots(idTypeNKnots);
end

% expType
optionsMFA.expType = [];
if isempty(optionsMFA.expType)
    optionsMFA.expType = expTypeAll{idExpType};
end

%最適化試行のID
optionsMFA.idOptimIter = [];
if isempty(optionsMFA.idOptimIter)
    optionsMFA.idOptimIter = idOptimIter;
end

%% データファイル名等の指定
preModelFileName = 'SimpleNetwork';
fileDirNames.expFileName = 'JamesData180205';
comptName = {'cell','media'};
fileDirNames.modelFileName = [preModelFileName];
fileDirNames.preSolFileName = ['preFinalResult' optionsMFA.expType num2str(optionsMFA.nKnots)];

% エクセルデータからモデルを更新するか？
% (urlreadが遺伝研で上手く使えない可能性があるため)
optionsMFA.isUpdateXlsNetwork = true;
% optionsMFA.isUpdateXlsNetwork = false;

%% popSizeと世代数
% popSizeについて、CMA-ESのデフォルトは(4 + floor(3*log(nParam)))
% nParam=30のとき、popSizeは14
% nParam = 50のとき、popSizeは15
% nParam=200のとき、popSizeは19
optionsOptimMH.popSize = 50;
optionsOptimMH.nMaxEval = 1*10^4;
% optionsOptimMH.nMaxEval = 0;

%% 出力フォルダ名
if ispc %Windows
    fileDirNames.saveDirName = [pwd '\ResultTest'];
else %Mac or UNIX
    fileDirNames.saveDirName = [pwd '/ResultNIG' optionsMFA.expType num2str(optionsMFA.nKnots)];
end

%% nKnotsのロード
% loadMFAInputでnKnotsの情報が必要になるのでここでロードする。
if optionsMFA.nKnots<0
    fileDirNames.bestNKnots = ['bestNKnots' optionsMFA.expType];
    loadData = load(fileDirNames.bestNKnots);
    
    optionsMFA.nKnotsRxns = loadData.bestNKnotsRxn;
    optionsMFA.nKnotsMets = loadData.bestNKnotsMet;
    optionsMFA.nKnotsNetRxns = loadData.bestNKnotsNetRxn;
    optionsMFA.nKnots = max(optionsMFA.nKnotsNetRxns);
    optionsMFA.isRxnSpecificNKnots = true;
    optionsMFA.patternNKnots = unique(optionsMFA.nKnotsNetRxns);
else
    optionsMFA.nKnotsRxns = []; % モデルを構築してから修正する
    optionsMFA.nKnotsMets = [];
    optionsMFA.nKnotsNetRxns = [];
    optionsMFA.isRxnSpecificNKnots = false;    
    optionsMFA.patternNKnots = optionsMFA.nKnots;
end

%% 培地情報
optionsMFA.isAccountMediaDrop = false;
% optionsMFA.isAccountMediaDrop = true;
mediaInfo.idCompt = 2;
mediaInfo.mets = {'Lac_ex'};
mediaInfo.initConcs = 10*10^3; % nM (= 10 uM ~ 1 uM)
% mediaInfo.initGlc = 10*10^6; % nM (= 10 mM)
% mediaInfo.initLac = 10^3; % nM (= 1 uM ~ 0 nM)
mediaInfo.mediaVol = 1.5*10^-3; % L
mediaInfo.cellAmount = 0.8; %mg
mediaInfo.mediaDropPerSampling = 0.08*10^-3; % mL, サンプリングに伴う培養液量の減少
optionsMFA.mediaInfo = mediaInfo;


%% 評価関数の指定
optionsMFA.objfunType = 1; % variance-weighted RSS of conc and MDV
% optionsMFA.objfunType = 2; % variance-weighted RSS of concMDV

%% compartmentの補正
% idComptCorr: コンパートメントを補正するのはどのコンパートメントか？
% コンパートメント1に対する比率をパラメータにする。
% optionsMFA.idComptCorrParam = [];
optionsMFA.idComptCorrParam = 2;
% optionsMFA.isParamCoefCorrCompt = true;
% optionsMFA.isParamCoefCorrCompt = false;

%% inner problemで解くcompartmentの指定
% optionsMFA.idComptInnerOptim = 1; % Lacをouter problemでフィット
optionsMFA.idComptInnerOptim = [1,2]; % Lacをinner problemでフィット

%% 指定するCVを使う場合
optionsMFA.isUseInputCV = false;
% optionsMFA.isUseInputCV = true;
optionsMFA.inputConcsCV = 0.1; %代謝物濃度と同位体濃度のCV
optionsMFA.inputMDVsCV = 0.1; %同位体割合のCV

%% フラックスの制限
optionsMFA.isFluxConstr = false;
% optionsMFA.isFluxConstr = true;
if   strncmp(fileDirNames.modelFileName,'G6P', 3)
    optionsMFA.isFluxConstr = false;
end
optionsMFA.isFluxConstr = false;
% optionsMFA.isFluxConstr = true;
optionsMFA.constr.fluxes.rxnNames = ...
    [repmat({'G6PD'}, 1, optionsMFA.nKnots+2), ...
    repmat({'PGD'}, 1, optionsMFA.nKnots+2),...
    repmat({'PRPS'}, 1, optionsMFA.nKnots+2),...
    repmat({'GPAT_ACL'}, 1, optionsMFA.nKnots+2)];
%     repmat({'GPAT'}, 1, optionsMFA.nKnots+2)];
optionsMFA.constr.fluxes.lb = zeros(1,(optionsMFA.nKnots+2)*4)+0.0001;
optionsMFA.constr.fluxes.ub = [...
    zeros(1,(optionsMFA.nKnots+2)*3)+2, ...
    zeros(1,(optionsMFA.nKnots+2)*1)+200];
optionsMFA.constr.fluxes.idFullKnot = repmat(1:optionsMFA.nKnots+2,1,4);

% optionsMFA.isFluxConstr = true;
% optionsMFA.constr.fluxes.rxnNames = repmat({'PRPS'}, 1, optionsMFA.nKnots+2);
% optionsMFA.constr.fluxes.lb = zeros(1,optionsMFA.nKnots+2)+0.01;
% optionsMFA.constr.fluxes.ub = zeros(1,optionsMFA.nKnots+2)+5;
% optionsMFA.constr.fluxes.idFullKnot = 1:optionsMFA.nKnots+2;

%% 初期代謝物濃度の制限
% 測定代謝物の初期代謝物濃度の制限については未対応 (knots2pCoefDataの修正が必要)
%GluとAspを評価する場合はそれらの初期濃度の制限も加える
% lbInitConcsについて、Ins, Glu:  mean([36,70,63]), Asp: mean([88,106,125]),
% Basは特になし

optionsMFA.isInitConcConstr = false;
% optionsMFA.isInitConcConstr = true;

if ~isempty(strfind(fileDirNames.modelFileName, 'Glycolysis'))
    optionsMFA.constr.initConcs.mets = {'Glc','Lac_ex'};
    optionsMFA.constr.initConcs.lb =  [1, 1];
    optionsMFA.constr.initConcs.ub =  [10^2, 10];
    
elseif ~isempty(strfind(fileDirNames.modelFileName, 'Central'))
%     optionsMFA.constr.initConcs.mets = {'Glc','Lac_ex', 'Mal_Fum', 'NAD'};
%     optionsMFA.constr.initConcs.lb =  [1, 1, 1, 1];
%     optionsMFA.constr.initConcs.ub =  [10^2, 10, 10^2, 10^2];
    optionsMFA.constr.initConcs.mets = {'Glc','Lac_ex', 'Mal_Fum'};
    optionsMFA.constr.initConcs.lb =  [1, 0.1, 1];
    optionsMFA.constr.initConcs.ub =  [10^2, 10^1, 10^2];
%     optionsMFA.constr.initConcs.mets = {'Glc','Lac_ex', 'Mal_Fum', 'AcCoA', 'AcCoA_cyto'};
%     optionsMFA.constr.initConcs.lb =  [1, 1, 1, 0.01, 0.01];
%     optionsMFA.constr.initConcs.ub =  [10^2, 10, 10^2, 1,1];
end

%% 代謝物濃度の変化速度の制限
optionsMFA.isConcRateConstr = false;
% optionsMFA.isConcRateConstr = true;
if   strncmp(fileDirNames.modelFileName,'G6P', 3)
    optionsMFA.isConcRateConstr = false;
end

if ~isempty(strfind(fileDirNames.modelFileName, 'Glycolysis'))
optionsMFA.constr.concRates.mets = ...
    {'Glc','G6P_F6P','Lac', 'Lac_ex', 'Glyc3P', 'Ala'};
optionsMFA.constr.concRates.lb =  -[10^2,10^2,10^2,10^2, 10^2, 10^2];
optionsMFA.constr.concRates.ub =  [10^2,10^2,10^2,10^2,10^2, 10^2];

elseif ~isempty(strfind(fileDirNames.modelFileName, 'Central'))
optionsMFA.constr.concRates.mets = ...
    {'Glc','G6P_F6P','Lac', 'Lac_ex', 'Mal_Fum', 'Glyc3P', 'Cit_Acon_IsoCit', 'Ala'};
optionsMFA.constr.concRates.lb =  -[10^2,10^2,10^2,10^2, 10^2, 10^2, 10^2, 10^2];
optionsMFA.constr.concRates.ub = [10^2,10^2,10^2,10^2,10^2, 10^2, 10^2, 10^2];
end

%% 独立フラックスの指定
optionsMFA.inputIndRxnNames = [];
% optionsMFA.inputIndRxnNames = {'G6PD'};

%% meta-heuristic最適化のオプション
optionsOptimMH.method = 'CMA-ES';
%         optionsGA.CMA.active = [];
optionsOptimMH.CMA.active = 1;
%         optionsGA.CMA.active = 2;
% OPTS.CMA.active = 1 turns on "active CMA" with a negative update
% of the covariance matrix and checks for positive definiteness.
% OPTS.CMA.active = 2 does not check for pos. def. and is numerically
% faster. Active CMA usually speeds up the adaptation and might
% become a default in near future.

optionsMFA.optionsOptimMH = optionsOptimMH;

%% local最適化のオプション
% local最適化の方法
% optionsMFA.optionsOptimLocal.method = 'lsqnonlin';  %テストしていない。
optionsMFA.optionsOptimLocal.method = 'fmincon';

%local最適化を実行するか？
optionsMFA.optionsOptimLocal.isPerform = true;
% optionsMFA.optionsOptimLocal.isPerform = false;

% 自作のヤコビ行列を使うか
optionsMFA.optionsOptimLocal.isUseMyJacobian= false;
% optionsMFA.optionsOptimLocal.isUseMyJacobian= true;

% local最適化の際にパラメータの範囲を広げる割合
optionsMFA.paramBoundExpandRatio = 100/99;
% optionsMFA.paramBoundExpandRatio = 2;

% local最適化での最小の評価回数
optionsMFA.optionsOptimLocal.nMaxEvalLocalOpt = 1*10^3;

%% パラメータの範囲のオプション
optionsMFA.lb.initConcs = 1*10^-1;
optionsMFA.lb.knotConcRates = -2*10^-1;% AspとかGluがないのであれば、小さくてよいと思う。
optionsMFA.lb.mesInitConcs = 1*10^-1;
optionsMFA.lb.mesKnotConcRates = -2*10^0;% AspとかGluがないのであれば、小さくてよいと思う。
optionsMFA.lb.knotFluxes = 10^-2;
optionsMFA.lb.knots =  0;
optionsMFA.lb.coefCorrCompt = 1;
% optionsMFA.lb.coefCorrCompt = 1-10^-3;
% optionsMFA.lb.concsOutCalib =  0.01;

optionsMFA.ub.initConcs = 1*10^2; % AspとかGluがないのであれば、小さくてよいと思う。
optionsMFA.ub.knotConcRates = +2*10^-1; % AspとかGluがないのであれば、小さくてよいと思う。
optionsMFA.ub.mesInitConcs = 1*10^2; % AspとかGluがないのであれば、小さくてよいと思う。
optionsMFA.ub.mesKnotConcRates = +2*10^0; % AspとかGluがないのであれば、小さくてよいと思う。
optionsMFA.ub.knotFluxes = 2*10^2;  % フラックスの最大値
% optionsMFA.ub.knotFluxes = 5*10^1;
% optionsMFA.ub.knots =  10; % knot timeの制限。（あまり好ましくはないが）
optionsMFA.ub.knots =  60;
optionsMFA.ub.coefCorrCompt = 1;
% optionsMFA.ub.coefCorrCompt = 1+10^-3;
% optionsMFA.ub.concsOutCalib =  2*10^5;


%% 関数名
optionsMFA.fxnName.selectIndVars = 'selectIndVars';
optionsMFA.fxnName.optimizationInMIMFA = 'optimizationInMIMFA';
optionsMFA.fxnName.prepQpModel = 'prepQpModel';
optionsMFA.fxnName.param2vector = 'param2vector';
optionsMFA.fxnName.prepConvertMat = 'prepConvertMat';
optionsMFA.fxnName.prepConvertMatDepKnots = 'prepConvertMatDepKnots';
optionsMFA.fxnName.createInitParams = 'createInitParams';
optionsMFA.fxnName.modifyQPModelObjFun = 'modifyQPModelObjFun';
optionsMFA.fxnName.smplParamsUnderLinConstr  = 'smplParamsUnderLinConstr';
optionsMFA.fxnName.solveMDVDynamics = 'solveMDVDynamics';
optionsMFA.fxnName.solveEMUDynamics = 'solveEMUDynamics';
optionsMFA.fxnName.vector2param = 'vector2param';
optionsMFA.fxnName.paramInd2LocalFull = 'paramInd2LocalFull';
optionsMFA.fxnName.saveMFASol = 'saveMFASol';

%% ODEのオプション
optionsODE.solver = 'ode15s';
% optionsODE.EMUTol= 10^-2;
optionsODE.AbsTol = 5*10^-4;
optionsODE.RelTol = 5*10^-5;
optionsODE.BDF = 'off';
optionsODE.MaxOrder = 3;
optionsODE.isUseMyJacobian = false;
optionsODE.isUseMyJPattern = true;
% optionsODE.Vectorized = 'off';

%odesetで指定するoption
optionsODE.activeOptionNames = {'AbsTol', 'RelTol', 'BDF', 'MaxOrder'};
optionsMFA.optionsODE = optionsODE;

%% ODEの検討
% ODEの検討を行うか？最適化初期パラメータをサンプリングしたのちに、他のODE条件を検討する。
% optionsMFA.isCheckODE= true;
optionsMFA.isCheckODE = false;

if optionsMFA.isCheckODE
    idODEOption = 1;
    while true
        idODEOption = idODEOption + 1;
        tmpOptionsODE = optionsODE;
        switch idODEOption
            case 2  % 以前のバージョン
                tmpOptionsODE.AbsTol = 1*10^-3;
                tmpOptionsODE.RelTol = 1*10^-4;
            case 3   %現在のベスト
                tmpOptionsODE.AbsTol = 5*10^-4;
                tmpOptionsODE.RelTol = 5*10^-5;
            otherwise
                break
        end
        optionsMFA.optionsODE(idODEOption) = tmpOptionsODE;
    end
end

%% その他の入力情報

% 以前の最適化結果を使うか？
% optionsMFA.isUsePreSol = true;
optionsMFA.isUsePreSol = false;

% 初期値 (0min) は定常にするか？
% optionsMFA.isInitSteadyState = true; %現在未対応。細胞外乳酸の初期値定常はイマイチ?,
optionsMFA.isInitSteadyState = false;

% 以前の最適化結果を使う場合のsubSolのID
optionsMFA.idPreSol = optionsMFA.idOptimIter;

% RSS計算に使うreplicate数
optionsMFA.nRepThreshold = 2;

% 評価するconcs, concMDVsのSDの閾値
optionsMFA.minConcsCV = 0;
% optionsMFA.minConcsCV = 0.01;

% 評価するMDVのSDの閾値
optionsMFA.minMDVsSD = 0;
% optionsMFA.minMDVsSD = 0.005;

% QP制約条件の時間step数
optionsMFA.nStepTimeConstrQP = 31;

% knot時間の最小の間隔
optionsMFA.minKnotTimeDiff = 1.01; % PWACentral287_2でminKnotTimeDiff = 1.1の場合を検討中。

% sumEMUの許容誤差
optionsMFA.tolEmuError= 10^-2;

% paramMHサンプリングとQPの優先順位(1: paramMHサンプリング, 2: QP)
optionsMFA.priorityParamMHSmplQP = 1;
% optionsMFA.priorityParamMHSmplQP = 2; % rxnSpecificNKnotsではこちらは使えない (concRatesとindFluxの間に制約条件があるため)

% paramMHにmesMetsのinitConcsを含めるか？
optionsMFA.isIncludeAllInitConcsInParamMH = true;
% optionsMFA.isIncludeAllInitConcsInParamMH = false;

% シミュレーション値に天然同位体補正を加えるか？
optionsMFA.isCorrectNatAbundanceSim = false;
% optionsMFA.isCorrectNatAbundanceSim = true;

% metaheuristic optimizationでQPを使うか？
optionsMFA.isUseQPInMH = true;
% optionsMFA.isUseQPInMH = false;

%% 入力データまとめ
inputData.fileDirNames = fileDirNames;
inputData.comptName= comptName;
inputData.optionsMFA = optionsMFA;

end

% %% meta-heuristic最適化のオプション
% % optionsMFA.GAModelType = 'SGA';
% % optionsGA.GAModelType = 'MGG';
% % optionsMFA.GAModelType = 'JGG';
% optionsOptimMH.method = 'CMA-ES';
% switch optionsOptimMH.method
%     case {'SGA'}
%         optionsOptimMH.selReproType = 'rouletteElite';
%         optionsOptimMH.selSurvType = 'allChildren';
%     case {'MGG'}
%         optionsOptimMH.selReproType = 'randomWOReplacement';
%         optionsOptimMH.selSurvType = 'rouletteElite';
% %         optionsMFA.nChildren = 2*optionsMFA.popSize;
%     case {'JGG'}
%         optionsOptimMH.selReproType = 'randomWOReplacement';
%         optionsOptimMH.selSurvType = 'eliteChildren';
% %         optionsMFA.nChildren = 2*optionsMFA.popSize;        
%     case {'CMA-ES'} 
% %         optionsGA.CMA.active = [];
%         optionsOptimMH.CMA.active = 1;
% %         optionsGA.CMA.active = 2;
%         % OPTS.CMA.active = 1 turns on "active CMA" with a negative update
%         % of the covariance matrix and checks for positive definiteness.
%         % OPTS.CMA.active = 2 does not check for pos. def. and is numerically
%         % faster. Active CMA usually speeds up the adaptation and might
%         % become a default in near future.
% 
%     otherwise
%         error([optionsOptimMH.GAModelType ' is no applicable'])
% end
% 
% switch optionsOptimMH.method
%     case {'SGA', 'MGG', 'JGG'}
%         % optionsGA.crossoverType = 'CMD-REX-U';
%         % optionsGA.crossoverType = 'REX-U';
%         % optionsGA.crossoverType = 'REX-N';
%         % optionsGA.crossoverType = 'UNDX-m';
%         optionsOptimMH.crossoverType = 'arithmetical';
%         optionsOptimMH.mutType = 'gaussian';
%         optionsOptimMH.mutCV = sqrt(1/10);
%         optionsOptimMH.mutRate = 0.1;
% end
% 
% optionsMFA.optionsOptimMH = optionsOptimMH;




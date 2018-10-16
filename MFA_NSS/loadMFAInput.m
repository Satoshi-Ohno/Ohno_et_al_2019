function inputData = loadMFAInput(arrayJobId)
%MFA�̓���

% 1����������̍œK���̎��s��
optionsMFA.nOptimIter = 10; % arrayJobId�ł��w��

% internal knots���̎�ސ�
typeNKnots = 0:2;
% typeNKnots = 2;
% typeNKnots = -1;  %�������Ƃ�nKnots���w�肷��
nTypeNKnots = length(typeNKnots);

%�����f�[�^�̎��
expTypeAll = {'Ins', 'Ctrl'};
nExpTypeAll = length(expTypeAll);

%% arrayJobID����internal knots�̐�, expType�̎w��
arrayJobIdsAll = 1:optionsMFA.nOptimIter*nTypeNKnots*nExpTypeAll;
matArrayJobIdsAll = reshape(arrayJobIdsAll, optionsMFA.nOptimIter, nTypeNKnots, nExpTypeAll);

for i = 1:nExpTypeAll
    if any(any(matArrayJobIdsAll(:,:,i)==arrayJobId))        
        idExpType = i;
        [idOptimIter, idTypeNKnots] = find(matArrayJobIdsAll(:,:,i)==arrayJobId);
        break
    end
end

% internal knots�̐��̎w��
optionsMFA.nKnots = []; 
% arrayJobId��nKnots���C��
if isempty(optionsMFA.nKnots)
    optionsMFA.nKnots = typeNKnots(idTypeNKnots);
end

% expType
optionsMFA.expType = [];
if isempty(optionsMFA.expType)
    optionsMFA.expType = expTypeAll{idExpType};
end

%�œK�����s��ID
optionsMFA.idOptimIter = [];
if isempty(optionsMFA.idOptimIter)
    optionsMFA.idOptimIter = idOptimIter;
end

%% �f�[�^�t�@�C�������̎w��
preModelFileName = 'SimpleNetwork';
fileDirNames.expFileName = 'JamesData180205';
comptName = {'cell','media'};
fileDirNames.modelFileName = [preModelFileName];
fileDirNames.preSolFileName = ['preFinalResult' optionsMFA.expType num2str(optionsMFA.nKnots)];

% �G�N�Z���f�[�^���烂�f�����X�V���邩�H
% (urlread����`���ŏ�肭�g���Ȃ��\�������邽��)
optionsMFA.isUpdateXlsNetwork = true;
% optionsMFA.isUpdateXlsNetwork = false;

%% popSize�Ɛ��㐔
% popSize�ɂ��āACMA-ES�̃f�t�H���g��(4 + floor(3*log(nParam)))
% nParam=30�̂Ƃ��ApopSize��14
% nParam = 50�̂Ƃ��ApopSize��15
% nParam=200�̂Ƃ��ApopSize��19
optionsOptimMH.popSize = 50;
optionsOptimMH.nMaxEval = 1*10^4;
% optionsOptimMH.nMaxEval = 0;

%% �o�̓t�H���_��
if ispc %Windows
    fileDirNames.saveDirName = [pwd '\ResultTest'];
else %Mac or UNIX
    fileDirNames.saveDirName = [pwd '/ResultNIG' optionsMFA.expType num2str(optionsMFA.nKnots)];
end

%% nKnots�̃��[�h
% loadMFAInput��nKnots�̏�񂪕K�v�ɂȂ�̂ł����Ń��[�h����B
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
    optionsMFA.nKnotsRxns = []; % ���f�����\�z���Ă���C������
    optionsMFA.nKnotsMets = [];
    optionsMFA.nKnotsNetRxns = [];
    optionsMFA.isRxnSpecificNKnots = false;    
    optionsMFA.patternNKnots = optionsMFA.nKnots;
end

%% �|�n���
optionsMFA.isAccountMediaDrop = false;
% optionsMFA.isAccountMediaDrop = true;
mediaInfo.idCompt = 2;
mediaInfo.mets = {'Lac_ex'};
mediaInfo.initConcs = 10*10^3; % nM (= 10 uM ~ 1 uM)
% mediaInfo.initGlc = 10*10^6; % nM (= 10 mM)
% mediaInfo.initLac = 10^3; % nM (= 1 uM ~ 0 nM)
mediaInfo.mediaVol = 1.5*10^-3; % L
mediaInfo.cellAmount = 0.8; %mg
mediaInfo.mediaDropPerSampling = 0.08*10^-3; % mL, �T���v�����O�ɔ����|�{�t�ʂ̌���
optionsMFA.mediaInfo = mediaInfo;


%% �]���֐��̎w��
optionsMFA.objfunType = 1; % variance-weighted RSS of conc and MDV
% optionsMFA.objfunType = 2; % variance-weighted RSS of concMDV

%% compartment�̕␳
% idComptCorr: �R���p�[�g�����g��␳����̂͂ǂ̃R���p�[�g�����g���H
% �R���p�[�g�����g1�ɑ΂���䗦���p�����[�^�ɂ���B
% optionsMFA.idComptCorrParam = [];
optionsMFA.idComptCorrParam = 2;
% optionsMFA.isParamCoefCorrCompt = true;
% optionsMFA.isParamCoefCorrCompt = false;

%% inner problem�ŉ���compartment�̎w��
% optionsMFA.idComptInnerOptim = 1; % Lac��outer problem�Ńt�B�b�g
optionsMFA.idComptInnerOptim = [1,2]; % Lac��inner problem�Ńt�B�b�g

%% �w�肷��CV���g���ꍇ
optionsMFA.isUseInputCV = false;
% optionsMFA.isUseInputCV = true;
optionsMFA.inputConcsCV = 0.1; %��ӕ��Z�x�Ɠ��ʑ̔Z�x��CV
optionsMFA.inputMDVsCV = 0.1; %���ʑ̊�����CV

%% �t���b�N�X�̐���
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

%% ������ӕ��Z�x�̐���
% �����ӕ��̏�����ӕ��Z�x�̐����ɂ��Ă͖��Ή� (knots2pCoefData�̏C�����K�v)
%Glu��Asp��]������ꍇ�͂����̏����Z�x�̐�����������
% lbInitConcs�ɂ��āAIns, Glu:  mean([36,70,63]), Asp: mean([88,106,125]),
% Bas�͓��ɂȂ�

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

%% ��ӕ��Z�x�̕ω����x�̐���
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

%% �Ɨ��t���b�N�X�̎w��
optionsMFA.inputIndRxnNames = [];
% optionsMFA.inputIndRxnNames = {'G6PD'};

%% meta-heuristic�œK���̃I�v�V����
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

%% local�œK���̃I�v�V����
% local�œK���̕��@
% optionsMFA.optionsOptimLocal.method = 'lsqnonlin';  %�e�X�g���Ă��Ȃ��B
optionsMFA.optionsOptimLocal.method = 'fmincon';

%local�œK�������s���邩�H
optionsMFA.optionsOptimLocal.isPerform = true;
% optionsMFA.optionsOptimLocal.isPerform = false;

% ����̃��R�r�s����g����
optionsMFA.optionsOptimLocal.isUseMyJacobian= false;
% optionsMFA.optionsOptimLocal.isUseMyJacobian= true;

% local�œK���̍ۂɃp�����[�^�͈̔͂��L���銄��
optionsMFA.paramBoundExpandRatio = 100/99;
% optionsMFA.paramBoundExpandRatio = 2;

% local�œK���ł̍ŏ��̕]����
optionsMFA.optionsOptimLocal.nMaxEvalLocalOpt = 1*10^3;

%% �p�����[�^�͈̔͂̃I�v�V����
optionsMFA.lb.initConcs = 1*10^-1;
optionsMFA.lb.knotConcRates = -2*10^-1;% Asp�Ƃ�Glu���Ȃ��̂ł���΁A�������Ă悢�Ǝv���B
optionsMFA.lb.mesInitConcs = 1*10^-1;
optionsMFA.lb.mesKnotConcRates = -2*10^0;% Asp�Ƃ�Glu���Ȃ��̂ł���΁A�������Ă悢�Ǝv���B
optionsMFA.lb.knotFluxes = 10^-2;
optionsMFA.lb.knots =  0;
optionsMFA.lb.coefCorrCompt = 1;
% optionsMFA.lb.coefCorrCompt = 1-10^-3;
% optionsMFA.lb.concsOutCalib =  0.01;

optionsMFA.ub.initConcs = 1*10^2; % Asp�Ƃ�Glu���Ȃ��̂ł���΁A�������Ă悢�Ǝv���B
optionsMFA.ub.knotConcRates = +2*10^-1; % Asp�Ƃ�Glu���Ȃ��̂ł���΁A�������Ă悢�Ǝv���B
optionsMFA.ub.mesInitConcs = 1*10^2; % Asp�Ƃ�Glu���Ȃ��̂ł���΁A�������Ă悢�Ǝv���B
optionsMFA.ub.mesKnotConcRates = +2*10^0; % Asp�Ƃ�Glu���Ȃ��̂ł���΁A�������Ă悢�Ǝv���B
optionsMFA.ub.knotFluxes = 2*10^2;  % �t���b�N�X�̍ő�l
% optionsMFA.ub.knotFluxes = 5*10^1;
% optionsMFA.ub.knots =  10; % knot time�̐����B�i���܂�D�܂����͂Ȃ����j
optionsMFA.ub.knots =  60;
optionsMFA.ub.coefCorrCompt = 1;
% optionsMFA.ub.coefCorrCompt = 1+10^-3;
% optionsMFA.ub.concsOutCalib =  2*10^5;


%% �֐���
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

%% ODE�̃I�v�V����
optionsODE.solver = 'ode15s';
% optionsODE.EMUTol= 10^-2;
optionsODE.AbsTol = 5*10^-4;
optionsODE.RelTol = 5*10^-5;
optionsODE.BDF = 'off';
optionsODE.MaxOrder = 3;
optionsODE.isUseMyJacobian = false;
optionsODE.isUseMyJPattern = true;
% optionsODE.Vectorized = 'off';

%odeset�Ŏw�肷��option
optionsODE.activeOptionNames = {'AbsTol', 'RelTol', 'BDF', 'MaxOrder'};
optionsMFA.optionsODE = optionsODE;

%% ODE�̌���
% ODE�̌������s�����H�œK�������p�����[�^���T���v�����O�����̂��ɁA����ODE��������������B
% optionsMFA.isCheckODE= true;
optionsMFA.isCheckODE = false;

if optionsMFA.isCheckODE
    idODEOption = 1;
    while true
        idODEOption = idODEOption + 1;
        tmpOptionsODE = optionsODE;
        switch idODEOption
            case 2  % �ȑO�̃o�[�W����
                tmpOptionsODE.AbsTol = 1*10^-3;
                tmpOptionsODE.RelTol = 1*10^-4;
            case 3   %���݂̃x�X�g
                tmpOptionsODE.AbsTol = 5*10^-4;
                tmpOptionsODE.RelTol = 5*10^-5;
            otherwise
                break
        end
        optionsMFA.optionsODE(idODEOption) = tmpOptionsODE;
    end
end

%% ���̑��̓��͏��

% �ȑO�̍œK�����ʂ��g�����H
% optionsMFA.isUsePreSol = true;
optionsMFA.isUsePreSol = false;

% �����l (0min) �͒��ɂ��邩�H
% optionsMFA.isInitSteadyState = true; %���ݖ��Ή��B�זE�O���_�̏����l���̓C�}�C�`?,
optionsMFA.isInitSteadyState = false;

% �ȑO�̍œK�����ʂ��g���ꍇ��subSol��ID
optionsMFA.idPreSol = optionsMFA.idOptimIter;

% RSS�v�Z�Ɏg��replicate��
optionsMFA.nRepThreshold = 2;

% �]������concs, concMDVs��SD��臒l
optionsMFA.minConcsCV = 0;
% optionsMFA.minConcsCV = 0.01;

% �]������MDV��SD��臒l
optionsMFA.minMDVsSD = 0;
% optionsMFA.minMDVsSD = 0.005;

% QP��������̎���step��
optionsMFA.nStepTimeConstrQP = 31;

% knot���Ԃ̍ŏ��̊Ԋu
optionsMFA.minKnotTimeDiff = 1.01; % PWACentral287_2��minKnotTimeDiff = 1.1�̏ꍇ���������B

% sumEMU�̋��e�덷
optionsMFA.tolEmuError= 10^-2;

% paramMH�T���v�����O��QP�̗D�揇��(1: paramMH�T���v�����O, 2: QP)
optionsMFA.priorityParamMHSmplQP = 1;
% optionsMFA.priorityParamMHSmplQP = 2; % rxnSpecificNKnots�ł͂�����͎g���Ȃ� (concRates��indFlux�̊Ԃɐ�����������邽��)

% paramMH��mesMets��initConcs���܂߂邩�H
optionsMFA.isIncludeAllInitConcsInParamMH = true;
% optionsMFA.isIncludeAllInitConcsInParamMH = false;

% �V�~�����[�V�����l�ɓV�R���ʑ̕␳�������邩�H
optionsMFA.isCorrectNatAbundanceSim = false;
% optionsMFA.isCorrectNatAbundanceSim = true;

% metaheuristic optimization��QP���g�����H
optionsMFA.isUseQPInMH = true;
% optionsMFA.isUseQPInMH = false;

%% ���̓f�[�^�܂Ƃ�
inputData.fileDirNames = fileDirNames;
inputData.comptName= comptName;
inputData.optionsMFA = optionsMFA;

end

% %% meta-heuristic�œK���̃I�v�V����
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




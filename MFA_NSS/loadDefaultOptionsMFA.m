function optionsMFA = loadDefaultOptionsMFA()

%% bounds
optionsMFA.lb.initConcs = 1*10^-1;
optionsMFA.lb.knotConcRates = -2*10^-1;
optionsMFA.lb.mesInitConcs = 1*10^-1;
optionsMFA.lb.mesKnotConcRates = -2*10^0;
optionsMFA.lb.knotFluxes = 10^-2;
optionsMFA.lb.knots =  0;
optionsMFA.lb.coefCorrCompt = 1;

optionsMFA.ub.initConcs = 1*10^2; 
optionsMFA.ub.knotConcRates = +2*10^-1; 
optionsMFA.ub.mesInitConcs = 1*10^2; 
optionsMFA.ub.mesKnotConcRates = +2*10^0;
optionsMFA.ub.knotFluxes = 2*10^2;
optionsMFA.ub.knots =  60;
optionsMFA.ub.coefCorrCompt = 1;

%% Compartment name
comptName = {'cell','media'};


%% nKnots
    optionsMFA.nKnotsRxns = [];
    optionsMFA.nKnotsMets = [];
    optionsMFA.nKnotsNetRxns = [];
    optionsMFA.isRxnSpecificNKnots = false;    

%% media information
optionsMFA.isAccountMediaDrop = false;

%% correction of compartment
% optionsMFA.idComptCorrParam = [];
optionsMFA.idComptCorrParam = 2;

%% set compartments which is solved by QP
optionsMFA.idComptInnerOptim = [1,2]; 

%% use input CV
optionsMFA.isUseInputCV = false;

%% constraint of flux
optionsMFA.isFluxConstr = false;

%% constraint of initial metabolite concentrations

optionsMFA.isInitConcConstr = false;

%% constraint of change in metabolite concentrations
optionsMFA.isConcRateConstr = false;

%% options for local optimizations
% local最適化の際にパラメータの範囲を広げる割合
optionsMFA.paramBoundExpandRatio = 100/99;
% optionsMFA.paramBoundExpandRatio = 2;

%% function name
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

%% options for ODEs
optionsODE.solver = 'ode15s';
optionsODE.AbsTol = 5*10^-4;
optionsODE.RelTol = 5*10^-5;
optionsODE.BDF = 'off';
optionsODE.MaxOrder = 3;

optionsODE.activeOptionNames = {'AbsTol', 'RelTol', 'BDF', 'MaxOrder'};
optionsMFA.optionsODE = optionsODE;

%% Others

% threshold of number of replicates for calculating RSS
optionsMFA.nRepThreshold = 2;

% minimum CV of metabolie concentrations
optionsMFA.minConcsCV = 0;

% minimum SD of MDV
optionsMFA.minMDVsSD = 0;
% optionsMFA.minMDVsSD = 0.005;

% number of times whe solving QP
optionsMFA.nStepTimeConstrQP = 31;

% mimimum difference of knot time
optionsMFA.minKnotTimeDiff = 1.01; % PWACentral287_2でminKnotTimeDiff = 1.1の場合を検討中。

% tolerance of EMU
optionsMFA.tolEmuError= 10^-2;

% Priolity
optionsMFA.priorityParamMHSmplQP = 1;
% optionsMFA.priorityParamMHSmplQP = 2; % rxnSpecificNKnotsではこちらは使えない (concRatesとindFluxの間に制約条件があるため)

% Whether initial metabolite concentrations are included in paramMH
optionsMFA.isIncludeAllInitConcsInParamMH = true;
% optionsMFA.isIncludeAllInitConcsInParamMH = false;

% Correct natural abundances
optionsMFA.isCorrectNatAbundanceSim = false;
% optionsMFA.isCorrectNatAbundanceSim = true;

% metaheuristic optimizationでQPを使うか？
optionsMFA.isUseQPInMH = true;
% optionsMFA.isUseQPInMH = false;

end
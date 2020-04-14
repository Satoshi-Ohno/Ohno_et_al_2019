function optionsMFA = loadDefaultOptionsMFA()

%% bounds
% optionsMFA.lb.concs = 1*10^-1;
% optionsMFA.lb.concRates = -2*10^-1;
% optionsMFA.lb.mesConcs = 1*10^-1;
% optionsMFA.lb.mesConcRates = -2*10^2;
% optionsMFA.lb.fluxes = 10^-2;
optionsMFA.lb.switchTimes =  0;

% optionsMFA.ub.concs = 1*10^1; 
% optionsMFA.ub.concRates = +2*10^-1; 
% optionsMFA.ub.mesConcs = 3*10^3; 
% optionsMFA.ub.mesConcRates = +2*10^2;
% optionsMFA.ub.fluxes = 1*10^2;
optionsMFA.ub.switchTimes =  60;

%% options for reaction dependent number of time intervals
% optionsMFA.isRxnDepNSwitchTimes = true;
optionsMFA.isRxnDepNSwitchTimes = false;
% optionsMFA.nSwitchTimesRxns = [];
% optionsMFA.nSwitchTimesMets = [];
% optionsMFA.nSwitchTimesNetRxns = [];
optionsMFA.nSwitchTimesRxns = [1;1;0;0;2;2;1;1;1;1];  % # of switch times for each reaction
optionsMFA.nSwitchTimesMets = [1;2;0;2;1]; % # of switch times for metabolites
optionsMFA.nSwitchTimesNetRxns = [1;1;0;0;2;2;1;1;1];% # of switch times for each net reaction

%% options for metaheuristic optimization
optionsOptimMH.isUseInitCovMat = true;

optionsMFA.optionsOptimMH = optionsOptimMH;

%% options for local optimizations
optionsMFA.paramBoundExpandRatio = 100/99;

%% function name
% optionsMFA.fxnName.selectIndVars = 'selectIndVars';
% optionsMFA.fxnName.optimizationInMIMFA = 'optimizationInMIMFA';
% optionsMFA.fxnName.prepQpModel = 'prepQpModel';
% optionsMFA.fxnName.param2vector = 'param2vector';
% optionsMFA.fxnName.prepConvertMat = 'prepConvertMat';
% optionsMFA.fxnName.prepConvertMatDepSwitchTimes = 'prepConvertMatDepSwitchTimes';
% optionsMFA.fxnName.createInitParams = 'createInitParams';
% optionsMFA.fxnName.modifyQPModelObjFun = 'modifyQPModelObjFun';
% optionsMFA.fxnName.solveMDVDynamics = 'solveMDVDynamics';
% optionsMFA.fxnName.solveEMUDynamics = 'solveEMUDynamics';
% optionsMFA.fxnName.vector2param = 'vector2param';
% optionsMFA.fxnName.paramInd2LocalFull = 'paramInd2LocalFull';

%% options for ODEs
optionsODE.solver = 'ode15s';
optionsODE.AbsTol = 5*10^-4;
optionsODE.RelTol = 5*10^-5;
optionsODE.BDF = 'off';
optionsODE.MaxOrder = 3;

optionsODE.activeOptionNames = {'AbsTol', 'RelTol', 'BDF', 'MaxOrder'};
optionsMFA.optionsODE = optionsODE;

%% Others
% number of times whe solving QP
optionsMFA.nStepTimeConstrQP = 31;

% mimimum difference of switch time
optionsMFA.minKnotTimeDiff = 1.01; 

% tolerance of EMU
optionsMFA.tolEmuError= 10^-2;

% Priolity
optionsMFA.priorityParamMHSmplQP = 1;

% Whether initial metabolite concentrations are included in paramMH
optionsMFA.isIncludeInitEndConcsInParamMH = true;
% optionsMFA.isIncludeInitEndConcsInParamMH = false;

% Correct natural abundances
optionsMFA.isCorrectNatAbundanceSim = false;
% optionsMFA.isCorrectNatAbundanceSim = true;

% use QP in metaheuristic optimization
optionsMFA.isUseQPInMH = true;
% optionsMFA.isUseQPInMH = false;

% Same initial flux and metabolite concentrations between Ins and Ctrl
% In the present version, metabolites in the Ins condition must be listed in the first half and metabolites in the Ctrl condition must be listed in the second half.
% The same is true for reactions.
% This will be modified in a future version so that reaction ID and metabolite ID can be used.
optionsMFA.isSameInitInsCtrl = false;

% Use metabolite concentrations at switch times as parameter during optimization
% optionsMFA.isUseConcAsParam = true;
optionsMFA.isUseConcAsParam = false;

% compartment IDs under steady state at initial time
optionsMFA.idComptInitSteadyState= [];  

% Use closed functions
% This is for a private setting. Closed functions are not open.
optionsMFA.isClosed = false;

end
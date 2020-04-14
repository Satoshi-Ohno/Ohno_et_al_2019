%% Example of metabolic flux analysis under non-steady-state for simple metabolic network

init;
import packFxnMFANSS.*

%% Options
% Excel file name of model
optionsMFA.xlsModelFileName = 'SimpleNetwork_InsCtrl';

% Excel file name of experimental data
optionsMFA.xlsExpDataFileName = '';  % NOT avaiable

% number of switch times
optionsMFA.nSwitchTimes = 1;

% initial population size for CMA-ES
optionsMFA.optionsOptimMH.initPopSize = 50;

% population size during CMA-ES
optionsMFA.optionsOptimMH.popSize = 5;
% optionsMFA.optionsOptimMH.popSize = 30;

% maximum number of evaluations in CMA-ES
optionsMFA.optionsOptimMH.nMaxEval = 5*10^1;
% optionsMFA.optionsOptimMH.nMaxEval =2*10^4;

% Negative update of covariance matrix
optionsMFA.optionsOptimMH.CMA.active = 1;

% mauximum number of evaluation in fmincon
optionsMFA.optionsOptimLocal.nMaxEvalLocalOpt = 5*10^1;
% optionsMFA.optionsOptimLocal.nMaxEvalLocalOpt = 1*10^4;

% Same initial flux and metabolite concentrations between Ins and Ctrl
% In the present version, metabolites in the Ins condition must be listed in the first half and metabolites in the Ctrl condition must be listed in the second half.
% The same is true for reactions.
% This will be modified in a future version so that reaction ID and metabolite ID can be used.
optionsMFA.isSameInitInsCtrl = true;  


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load model defined in Excel file
model.name = optionsMFA.xlsModelFileName;
model = loadXlsNetwork(model, optionsMFA);

%% Experimental data
load('expData_SimpleNetwork_InsCtrl')

%% Run metabolic flux analysis under non-steady-state
[solMFA, model, expData, optionsMFA] = run_MFA_NSS(model, expData, optionsMFA);

%% Draw estimated metabolite concentrations and mass isotopomer fractions
optionsPlot.plotMet = model.mets(optionsMFA.varSet.idNonPoolMets)';
figH = drawFigFittingResult(model, expData, solMFA, optionsMFA, optionsPlot);

%% Draw estimated fluxes
figH = drawFigFluxTimeCourse(model, expData, solMFA, optionsMFA, optionsPlot);


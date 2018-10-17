%% Example of metabolic flux analysis under non-steady-state for simple metabolic network

clc
clear 
close all
import packFxnMIMFA.*
import fxnDrawFigs.*

%% Excel file name for model

xlsModelFileName = 'SimpleNetwork';

%% Experimental data
load('expDataSimpleNetwork')

%% Options

% number of switch times
optionsMFA.nSwitchTimes = 1;

% population size in CMA-ES
optionsMFA.optionsOptimMH.popSize = 50;

% maximum number of evaluations in CMA-ES
optionsMFA.optionsOptimMH.nMaxEval = 1*10^4;

% Negative update of covariance matrix
optionsMFA.optionsOptimMH.CMA.active = 1;

% mauximum number of evaluation in fmincon
optionsMFA.optionsOptimLocal.nMaxEvalLocalOpt = 1*10^4;


%% Run metabolic flux analysis under non-steady-state
[solMFA, model, expData, optionsMFA] = run_MFA_NSS(xlsModelFileName, expData, optionsMFA);% keyboard;

%% Draw estimated metabolite concentrations and mass isotopomer fractions
figH = drawFigFittingResult(model, expData, solMFA, optionsMFA, options);

%% Draw estimated fluxes
figH = drawFigFluxTimeCourse(model, expData, solMFA, optionsMFA, options);


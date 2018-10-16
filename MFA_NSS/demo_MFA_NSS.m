%% Example of metabolic flux analysis under non-steady-state for simple metabolic network

clc
clear 
close all
import packFxnMIMFA.*
import fxnDrawFigs.*

%% Excel file name for model

xlsModelFileName = 'SimpleNetwork';
% load('modelSimpleNetwork')

%% Experimental data
load('expDataSimpleNetwork')

%% Options
optionsMFA.nKnots = 1;

% population size in CMA-ES
optionsOptimMH.popSize = 30;

% maximum number of evaluations in CMA-ES
optionsOptimMH.nMaxEval = 1*10^2;

% Negative update of covariance matrix
optionsOptimMH.CMA.active = 1;

optionsMFA.optionsOptimMH = optionsOptimMH;

% mauximum number of evaluation in fmincon
optionsMFA.optionsOptimLocal.nMaxEvalLocalOpt = 1*10^2;

%% Run metabolic flux analysis under non-steady-state

[solMFA, model, expData, optionsMFA] = run_MFA_NSS(xlsModelFileName, expData, optionsMFA);% keyboard;

%% Draw estimated metabolite concentrations and mass isotopomer fractions
options.plotMet = model.mets(optionsMFA(1).varSet.idNonPoolMets)';
options.typeFigStr = 'S2';
options.alphaCI = 0.9;
[tmpFigH] = drawFigFittingResult(model, expData, solMFA, [], optionsMFA, options);


%% Draw estimated fluxes
options.nPlotSol = min(1);
options.isRecalcErrorbar = false;
options.typePlotErrorbar = 2;
options.typeDrawCIArea = 2;
options.isPlotAllFitting = true;
options.typeFluxAll = 2;

legendList = {'Ins'};
plotColor = cell(1,length(options.typeFluxAll ));
plotColor{1} = 'rb';

options.legendList = legendList;
options.plotColor = plotColor;
[tmpFigH] = drawFigFluxTimeCourse(model, expData, solMFA, [], optionsMFA, options);


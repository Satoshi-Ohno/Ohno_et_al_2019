%% Example of identification of effective flux regulators

clc
clear 
close all

%% Flux and potential flux regulator data for Gys
% A part of phosphorylation and allosteric effectors are considered in this example

for i = 1 : 6
    switch i
        case 1
            expDataKA(i).regulatorType = 'flux';
            expDataKA(i).regulatorName = [];
            expDataKA(i).ScoefRegulatorMet = [];
            expDataKA(i).level = [0.8713,0.8406,0.8024,0.7474,0.5280,0.1504,0.1440,0.1360,0.1200,0.0559];
        case 2
            expDataKA(i).regulatorType = 'phospho';
            expDataKA(i).regulatorName = 'Gys_S728';
            expDataKA(i).ScoefRegulatorMet = [];
            expDataKA(i).level = [1.3200,3.5948,8.2426,12.1022,10.8810,1,1,1,1,1];
        case 3
            expDataKA(i).regulatorType = 'phospho';
            expDataKA(i).regulatorName = 'Gys_S412';
            expDataKA(i).ScoefRegulatorMet = [];
            expDataKA(i).level =[1.0912,1.1381,1.0882,1.3078,1.3925,1,1,1,1,1];
        case 4
            expDataKA(i).regulatorType = 'alloA';
            expDataKA(i).regulatorName = 'G6P';
            expDataKA(i).ScoefRegulatorMet = [];
            expDataKA(i).level =[2.5641,3.1500,4.3369,4.1904,4.9754,1.8616,1.1573,1.0714,1.0790,1.1788];
        case 5
            expDataKA(i).regulatorType = 'alloI';
            expDataKA(i).regulatorName = 'UTP';
            expDataKA(i).ScoefRegulatorMet = [];
            expDataKA(i).level =[7.3889,7.4198,7.8225,7.4960,7.7280,6.9506,7.1156,6.9694,6.3874,6.6384];
        case 6
            expDataKA(i).regulatorType = 'substrate';
            expDataKA(i).regulatorName = 'UDPGlc';
            expDataKA(i).ScoefRegulatorMet = 1;
            expDataKA(i).level =[5.6277,4.6723,3.9777,3.1110,2.7652,5.4453,5.4299,5.4126,4.8831,4.9040];
    end
end

%% Options

% time at which flux is estimated
optionsIEFR.timeKinetic = [1,5,10,20,60]; 

% Parameter is estimated using only net flux or both forward and backward fluxes
optionsIEFR.isEstForBack = false;  % only net flux
% optionsIEFR.isEstForBack = true;  % both forward and backward fluxes

%% Run identification of effective flux regulators
[solMinAIC, solNoReg, optionsIEFR] = run_IEFR(expDataKA, optionsIEFR);



%% Calculate regulation coefficient
[regCoef, regCoefFxn] = calcRegCoef(solMinAIC, optionsIEFR);

%% Draw input flux vs estimated flux
figH = drawScatterInputVsEstimatedFlux(solMinAIC, solNoReg, optionsIEFR);

%% Draw regulation coefficient
figH = drawRegCoef(regCoefFxn, optionsIEFR);

%% Draw time-averaged regulation coefficient
figH = drawTimeAveragedRegCoef(regCoefFxn, optionsIEFR);



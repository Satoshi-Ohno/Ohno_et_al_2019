%% Example of identification of effective flux regulators

clc
clear 
close all

%% Flux and potential flux regulator data for Pfk1
% A part of phosphorylation and allosteric effectors are considered in this example

for i = 1 : 7
    switch i
        case 1
            expDataKM(i).regulatorType = 'flux';
            expDataKM(i).regulatorName = [];
            expDataKM(i).ScoefRegulatorMet = [];
            expDataKM(i).level = [20.0395990337418,9.94715895287550,11.8658551409600,15.7032475171290,31.0528170217970,17.7953128925431,2.55485671990928,2.79628879491361,3.27915294492228,5.21060954495773];
        case 2
            expDataKM(i).regulatorType = 'phospho';
            expDataKM(i).regulatorName = 'Pfka_S667';
            expDataKM(i).ScoefRegulatorMet = [];
            expDataKM(i).level = [1.20085001615290,0.990430012329094,1.17955998677490,1.01103733329598,1.04898765347995,1,1,1,1,1];
        case 3
            expDataKM(i).regulatorType = 'alloA';
            expDataKM(i).regulatorName = 'F6P';
            expDataKM(i).ScoefRegulatorMet = [];
            expDataKM(i).level =[0.597228666666667,0.710745000000000,1.05396033333333,1.03268733333333,1.07962500000000,0.439317666666667,0.280586333333333,0.237791000000000,0.272813000000000,0.284649666666667];
        case 4
            expDataKM(i).regulatorType = 'alloA';
            expDataKM(i).regulatorName = 'F2,6P';
            expDataKM(i).ScoefRegulatorMet = [];
            expDataKM(i).level =[0.0719166477630962,0.0761529386723917,0.0797743359793533,0.0723492040443248,0.0535158023322439,0.0791670193963261,0.0805244591730452,0.0758469145964414,0.0753808187990119,0.0654039495521192];
        case 5
            expDataKM(i).regulatorType = 'alloI';
            expDataKM(i).regulatorName = 'ADP';
            expDataKM(i).ScoefRegulatorMet = [];
            expDataKM(i).level =[3.59604766666667,3.20079666666667,3.93733966666667,3.91533466666667,2.83915266666667,3.29844866666667,3.78804266666667,3.61985966666667,3.97364700000000,3.60177933333333];
        case 6
            expDataKM(i).regulatorType = 'substrate';
            expDataKM(i).regulatorName = 'F6P';
            expDataKM(i).ScoefRegulatorMet = 1;
            expDataKM(i).level =[0.597228666666667,0.710745000000000,1.05396033333333,1.03268733333333,1.07962500000000,0.439317666666667,0.280586333333333,0.237791000000000,0.272813000000000,0.284649666666667];
        case 7
            expDataKM(i).regulatorType = 'substrate';
            expDataKM(i).regulatorName = 'ATP';
            expDataKM(i).ScoefRegulatorMet = 1;
            expDataKM(i).level =[56.3593006666667,52.1779026666667,49.4896520000000,48.2689143333333,51.5076256666667,54.5921446666667,52.3016713333333,50.5915110000000,52.8876430000000,51.1218803333333];
    end
end

%% Options

% time at which flux is estimated
optionsKM.timeKinetic = [1,5,10,20,60]; 

% Parameter is estimated using only net flux or both forward and backward fluxes
optionsKM.isEstForBack = false;  % only net flux
% optionsIEFR.isEstForBack = true;  % both forward and backward fluxes

%% Run identification of effective flux regulators
[solMinAIC, solNoReg, optionsKM] = run_KM(expDataKM, optionsKM);



%% Calculate regulation coefficient
[regCoef, regCoefFxn] = calcRegCoef(solMinAIC, optionsKM);

%% Draw input flux vs estimated flux
figH = drawScatterInputVsEstimatedFlux(solMinAIC, solNoReg, optionsKM);

%% Draw regulation coefficient
figH = drawRegCoef(regCoefFxn, optionsKM);

%% Draw time-averaged regulation coefficient
figH = drawTimeAveragedRegCoef(regCoefFxn, optionsKM);



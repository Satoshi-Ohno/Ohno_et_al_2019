function optionsIEFR = loadDefaultOptionsIEFR()

%% Optimization

% iteration of optimization
optionsIEFR.nIterOpt = 2;

% parameter bounds
optionsIEFR.optionsOpt.maxAbsVal = 10^2;
optionsIEFR.optionsOpt.minAbsVal = 10^-6;

% population size in CMA-ES
optionsIEFR.optionsOpt.popSize = 10;

% maximum number of evaluations in CMA-ES
optionsIEFR.optionsOpt.nMaxEval = 1*10^3;

% score (RSS) when optimization problem is infreasible
optionsIEFR.optionsOpt.infeasibleScore = 10^7;



end

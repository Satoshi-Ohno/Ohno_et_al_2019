%% Cacluate score (RSS)
function [score, scoreMet, residual, resVec2MDV]= calcScore(model, expData, optionsMFA, sol)

[residual, chisqValMDV,resVec2MDV]=calcResidual(model, expData, optionsMFA, sol);
score = sum(residual.^2);
scoreMet= chisqValMDV;

end
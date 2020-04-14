%% Calculation of RSS
function [RSS, residuals, fluxEst, term] = calcRSS(param, fxnCalcRes, infeasibleScore)

[residuals, fluxEst, term] = fxnCalcRes(param);

tmpRes = residuals;
tmpRes(isinf(tmpRes))=0;
tmpRes(isnan(tmpRes))=0;
RSS = sum(tmpRes.^2);

if nnz(tmpRes)< length(tmpRes)/2
    RSS=infeasibleScore;
end
end



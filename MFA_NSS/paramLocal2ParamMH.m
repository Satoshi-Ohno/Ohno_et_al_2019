function paramMH = paramLocal2ParamMH(model, expData, optionsMFA, paramLocal)

transformMat = optionsMFA.transformMat;

%% paramLocal -> paramMH
paramMH = transformMat.paramLocal2ParamMH * paramLocal;

%% log transformation
paramMH(transformMat.isLogParamMH) = log10(paramMH(transformMat.isLogParamMH));

end
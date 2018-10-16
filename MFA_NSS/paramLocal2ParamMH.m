function paramMH = paramLocal2ParamMH(model, expData, optionsMFA, paramLocal)

field2var(optionsMFA.varSet)
convertMat = optionsMFA.convertMat;

%% paramLocal -> paramMH
paramMH = convertMat.paramLocal2ParamMH * paramLocal;

%% log transformation
paramMH(convertMat.isLogParamMH) = log10(paramMH(convertMat.isLogParamMH));

end
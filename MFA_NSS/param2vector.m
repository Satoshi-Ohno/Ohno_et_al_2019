%% vectorize parameters
function     paramVector = param2vector(...
    paramStruct, optionsMFA, optimType, isLogTransform, isIndParam)

if nargin <= 4
    isIndParam = false;
end
if nargin <= 3
    switch optimType
        case {'metaheuristic'}
            isLogTransform = true;
        case {'local', 'init'}
            isLogTransform = false;
    end
end
if nargin<=2
    optimType = 'metaheuristic';
    isLogTransform = true;
end

if isIndParam
    switch optimType
        case {'metaheuristic'}
            isIndParamsVec = optionsMFA.isIndParams.MH;
        case {'local', 'init'}
            isIndParamsVec = optionsMFA.isIndParams.local;
    end
end

if optionsMFA.isUseConcAsParam
    fieldNamesParam = {'switchTimeConcs', 'initConcRates', 'switchTimeFluxes', 'switchTimes'} ;
else
    fieldNamesParam = {'initConcs', 'switchTimeConcRates', 'switchTimeFluxes', 'switchTimes'} ;
end

%% vectorize parameters
for f = 1 : length(fieldNamesParam)
    tmpParam = paramStruct.(fieldNamesParam{f});
    vecParamStruct.(fieldNamesParam{f}) = reshape(...
        tmpParam, numel(tmpParam), 1);
end

%% log transformation
if isLogTransform
    for f = 1 : length(fieldNamesParam)
        if ~strcmp(fieldNamesParam(f), {'switchTimeConcRates', 'initConcRates'})
            vecParamStruct.(fieldNamesParam{f}) = log10(vecParamStruct.(fieldNamesParam{f}));
        end
    end
end

%% Merge parameters
paramVector = zeros(optionsMFA.idParamLocal.nParam,1);
ff =0;
for f = 1 : length(fieldNamesParam)
    nTmpParam = length(vecParamStruct.(fieldNamesParam{f}));
    paramVector(ff+(1:nTmpParam)) = vecParamStruct.(fieldNamesParam{f});
    ff = ff + nTmpParam;
end
paramVector = paramVector(1:ff);

%% Select independent parameters
if isIndParam
    paramVector = paramVector(isIndParamsVec);
end


end
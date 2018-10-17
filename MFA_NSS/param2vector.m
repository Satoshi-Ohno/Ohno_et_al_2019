%% vectorize parameters
function     paramVector = param2vector171206(...
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
            isIndParams = optionsMFA.isIndParams.MH;
        case {'local', 'init'}
            isIndParams = optionsMFA.isIndParams.local;
    end
end
fieldNamesParam = fieldnames(paramStruct);

%% vectorize parameters
for f = 1 : length(fieldNamesParam)
    tmpParam = paramStruct.(fieldNamesParam{f});
    vecParamStruct.(fieldNamesParam{f}) = reshape(...
        tmpParam, numel(tmpParam), 1);
end

%% log transformation
if isLogTransform
    for f = 1 : length(fieldNamesParam)
        if ~strcmp(fieldNamesParam(f), {'switchTimeConcRates'})
            vecParamStruct.(fieldNamesParam{f}) = log10(vecParamStruct.(fieldNamesParam{f}));
        end
    end
end

%% Merge parameters
paramVector = [...
    vecParamStruct.initConcs;...
    vecParamStruct.switchTimeConcRates;...
    vecParamStruct.switchTimeFluxes;...
    vecParamStruct.switchTimes;...
    ];

%% Select independent parameters
if isIndParam
    paramVector = paramVector(isIndParams);
end


end
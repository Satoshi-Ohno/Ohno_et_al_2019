%% Calculation of regulation coefficient
function [regCoef, regCoefFxn] = calcRegCoef(solMinAIC, optionsKM);


timeKinetic = optionsKM.timeKinetic;
nTimeKinetic = length(timeKinetic);

fxnRegulatorNames = {'phospho', 'alloAI', 'subsProdMet', 'unaccounted'};
legendList = {'\rho_{P}', '\rho_{A}', ...
    '\rho_{S}', '\rho_{U}'};

idList = 1:nTimeKinetic*2;
timeRho = timeKinetic(1:nTimeKinetic);

paramInfoAll = optionsKM.paramInfo;
paramInfo = solMinAIC.paramInfo;

%% Prepare output
regCoef.fxn = [];
regCoef.regulatorTypes = [];
regCoef.regulatorNames = [];
regCoef.fxnVal = [];
regCoef.rho = [];
regCoef.timeRho = [];

regCoefFxn.fxnRegulatorNames = fxnRegulatorNames;
regCoefFxn.rhoMat = zeros(length(fxnRegulatorNames), length(timeRho));
regCoefFxn.rhoMat(end,:) = 1;
regCoefFxn.rhoMatTimeAverage = zeros(length(fxnRegulatorNames), 1);
regCoefFxn.rhoMatTimeAverage(end) = 1;
regCoefFxn.rhoMR = zeros(1, length(timeRho));
regCoefFxn.timeRho =timeRho;
regCoefFxn.legendList = legendList;


%% Calcuation of each term
residual = solMinAIC.residuals(idList);
tmpFieldNames = fieldnames(solMinAIC.factor);
tmpFieldNames = [{'k1'}; setdiff(tmpFieldNames, {'fS_b', 'k1'})];
for f = 1 : length(tmpFieldNames)
    switch tmpFieldNames{f}
        case {'k1'}
            fxnName ='fU';
            regulatorTypes = {'unaccounted'};
            regulatorNames= {'unaccounted'};
            fxnVal = exp(residual');
        case {'fS_f'}
            fxnName = 'fS';
            regulatorTypes = {'subsProdMet'};
            regulatorNames = {'subsProdMet'};
            if isempty(solMinAIC.factor.fS_b)
                fxnVal = solMinAIC.factor.fS_f(:, idList);
            else
                fxnVal = solMinAIC.factor.fS_f(:, idList)-solMinAIC.factor.fS_b(:, idList);
            end
        otherwise
            fxnName = tmpFieldNames{f};
            isRegParam = solMinAIC.paramInfo.isRegParam;
            regulatorTypes = paramInfoAll.paramRegulatorTypes(paramInfo.paramRegulatorTypeIds(isRegParam));
            regulatorNames = paramInfoAll.paramRegulatorNames(paramInfo.paramRegulatorNameIds(isRegParam));
            fxnVal = solMinAIC.factor.(fxnName)(:, idList);
    end
    regCoef(f).fxn = fxnName;
    regCoef(f).regulatorTypes = regulatorTypes;
    regCoef(f).regulatorNames = regulatorNames;
    regCoef(f).fxnVal = fxnVal;
    regCoef(f).rho = [];
    regCoef(f).timeRho = timeRho;
    regCoef(f).deltaLogFxnVal = [];
end

%% Calculate regulation coefficient
for f =  1:length(regCoef)
    regCoef(f).deltaLogFxnVal = log(abs(regCoef(f).fxnVal(:,1:nTimeKinetic)))...
        - log(abs(regCoef(f).fxnVal(:,nTimeKinetic+(1:nTimeKinetic))));
end

deltaLogFxnVal = vertcat(regCoef.deltaLogFxnVal);
tmpSquareSum = sum(deltaLogFxnVal.^2,1);
for f =  1:length(regCoef)
    regCoef(f).rho =  regCoef(f).deltaLogFxnVal.^2 * diag(1./tmpSquareSum);
end
loc =   strcmp({regCoef.fxn}, {'fU'});
regCoef(loc).rho = 1-sum(vertcat(regCoef(~loc).rho),1);
regCoef = rmfield(regCoef, 'deltaLogFxnVal');

%% Sum regulation coefficient among the same regulator types
tmpRegCoefMat = vertcat(regCoef.rho);
regCoefRegulatorTypes = vertcat(regCoef.regulatorTypes);
fxnRhoMat = zeros(length(fxnRegulatorNames), size(tmpRegCoefMat,2));
for f = 1 : length(fxnRegulatorNames)
    switch fxnRegulatorNames{f}
        case {'phospho'}
            loc = strcmp(regCoefRegulatorTypes , 'phospho');
        case {'alloAI'}
            loc = strncmp(regCoefRegulatorTypes , 'allo', 4);
        case {'subsProdMet'}
            loc = strcmp(regCoefRegulatorTypes , 'subsProdMet');
        case {'unaccounted'}
            loc = strcmp(regCoefRegulatorTypes , 'unaccounted');
    end
    fxnRhoMat(f,:) = sum(tmpRegCoefMat(loc,:),1);
end

rhoMatTimeAverage = zeros(length(fxnRegulatorNames),1);
for f = 1 : length(fxnRegulatorNames)
    for t = 1 : nTimeKinetic -1
        rhoMatTimeAverage(f) = rhoMatTimeAverage(f)+(fxnRhoMat(f,t)+fxnRhoMat(f,t+1))*diff(timeKinetic(t:t+1))/2;
    end
    rhoMatTimeAverage(f) = rhoMatTimeAverage(f)/(timeKinetic(end)-timeKinetic(1));
end

regCoefFxn.fxnRegulatorNames = fxnRegulatorNames;
regCoefFxn.rhoMat = single(fxnRhoMat);
regCoefFxn.rhoMatTimeAverage = single(rhoMatTimeAverage);
regCoefFxn.timeRho =regCoef(1).timeRho;
regCoefFxn.legendList = legendList;

end
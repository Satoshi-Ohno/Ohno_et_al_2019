%% Calculation of residuals
function [residual, scoreMet, resVec2MDV]= calcResidual(model, expData, optionsMFA, sol)

field2var(optionsMFA.varSet)

idCompTimeSim = zeros(1,length(expData.time));
for t = 1 : nExpTime
    [~,idCompTimeSim(t)] = min(abs(sol.timeSim-expData.time(t)));
end

nMaxVec = nEvalMDVMets*nExpTime*(max(model.metInfo.nCarbonMets(idEvalMDVMets))+1) ... % MDV
    + nEvalConcMets*nExpTime ;% Conc

%% Make vector of MDV or concMDV

vecExp = zeros(nMaxVec,1);
vecSim = zeros(nMaxVec,1);
vecSD = zeros(nMaxVec,1);
scoreMet = [];
resVec2MDV.colLabel = {'idMet', 'idTime', 'idMDV'};
resVec2MDV.id = zeros(nMaxVec,1);
mm = 0;
for m = 1 : nEvalMDVMets
    idMet = idEvalMDVMets(m);
    nCarbon = model.metInfo.nCarbonMets(idMet);
    isUseMass = model.metInfo.useMassMat(idMet,1:nCarbon+1);
    
    isComp = expData.MDVs{idMet}>0;  
    isComp(~isUseMass, :) = 0;
    nComp = (nCarbon+1) * nExpTime;
    
    % Vector of measuremets
    vecExp(mm+(1:nComp)) = reshape(expData.MDVs{idMet}, nComp, 1);
    
    % Normalization by the sum of MDV to be compared
    tmpMDVsSim =  sol.MDVsSim{idMet}(:,idCompTimeSim);
    tmpMDVsSim(~isComp) = 0;
    sumCompMDVsSim = sum(tmpMDVsSim,1);
    tmpMDVsExp = expData.MDVs{idMet};
    tmpMDVsExp(~isComp) = 0;
    sumCompMDVsExp = sum(tmpMDVsExp,1);
    tmpMDVsSim = tmpMDVsSim ./ repmat(sumCompMDVsSim./sumCompMDVsExp, nCarbon+1,1) ;
    
    % Vector of estimates
    vecSim(mm+(1:nComp)) = reshape(tmpMDVsSim, nComp, 1);
    
    % Vector of standard deviations
    vecSD(mm+(1:nComp)) = reshape(expData.MDVsSD{idMet}, nComp, 1);
    
    % set SD to nan for data not to be compared
    tmpVecSD = vecSD(mm+(1:nComp));
    isCompVec = reshape(isComp, nComp, 1);
    tmpVecSD(~isCompVec) = nan;
    vecSD(mm+(1:nComp)) = tmpVecSD;
        
    % index of residulas
    resVec2MDV.id(mm+(1:nComp),1) = idMet;  % met
    [idMDV, idTime] = find(true(nCarbon+1, nExpTime));
    resVec2MDV.id(mm+(1:nComp),2) = idTime;  % time
    resVec2MDV.id(mm+(1:nComp),3) = idMDV;  % MDV
    
    mm = mm + nComp;
end

%% Make vector of metabolite concentrations
for m = 1 : nEvalConcMets
    idMet = idEvalConcMets(m);
    
    %identify time to be compared
    isComp = expData.concs(idMet,:)>0;  
    nComp = nExpTime;
    
    % vector of measurements
    vecExp(mm+(1:nComp)) = reshape(expData.concs(idMet,:), nComp, 1);
    
    % vector of estimates
    tmpConcsSim = sol.concs(idMet,idCompTimeSim);
    
    % correct according to coefCorrCompt
    idCompt = model.metInfo.compt(idMet);
    coefCorrCompt = [1;sol.coefCorrCompt];
    tmpTimeSim = sol.timeSim;
    if isempty(sol.coefCorrCompt)
        if idCompt>=2
            tmpConcsSim = corrComptCell2Media(model, expData, optionsMFA, sol.concs(idMet,:), tmpTimeSim);
        else
            tmpConcsSim = tmpConcsSim* coefCorrCompt(idCompt);
        end
    else
        tmpConcsSim = tmpConcsSim* coefCorrCompt(idCompt);
    end    
    vecSim(mm+(1:nComp)) = reshape(tmpConcsSim, nComp, 1);
    
    % vector of standard deviations
    vecSD(mm+(1:nComp)) = reshape(expData.concsSD(idMet,:), nComp, 1);
    
    % set SD to nan for data not to be compared
    tmpVecSD = vecSD(mm+(1:nComp));
    isCompVec = isComp;
    tmpVecSD(~isCompVec) = nan;
    vecSD(mm+(1:nComp)) = tmpVecSD;
    
    % index of residulas
    resVec2MDV.id(mm+(1:nComp),1) = idMet;  % met
    idTime = 1:nExpTime;
    resVec2MDV.id(mm+(1:nComp),2) = idTime;  % time
    resVec2MDV.id(mm+(1:nComp),3) = 0;  % MDV
    
    mm = mm + nComp;
end
vecExp = vecExp(1:mm);
vecSim = vecSim(1:mm);
vecSD = vecSD(1:mm);
resVec2MDV.id = resVec2MDV.id(1:mm,:);

isEval = ~isnan(vecSD);
vecExp = vecExp(isEval);
vecSim = vecSim(isEval);
vecSD = vecSD(isEval);
resVec2MDV.id = resVec2MDV.id(isEval,:);


%% Calculate residuals
residual = (vecExp-vecSim)./vecSD;

end

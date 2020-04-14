%% Make vector of measured metabolite concentrations
function vecExpConcStruct = makeVecExpConcs(model, expData, optionsMFA)

field2var(optionsMFA.varSet);

nExpConcs  = nEvalConcMets*nExpTime;
vecExpConcs = zeros(nExpConcs,1);
vecSDs = zeros(nExpConcs,1);
vecIdMet = zeros(nExpConcs,1);
vecIdTime = zeros(nExpConcs,1);
vecIsComp = zeros(nExpConcs,1);
mm = 0;
for m = 1 : nEvalConcMets
    idMet = idEvalConcMets(m);
        
    isComp = expData.concs(idMet,:)>0;  
    
    tmpExpConcs = reshape(expData.concs(idMet, :), nExpTime, 1);
    tmpSDs = reshape(expData.concsSD(idMet, :), nExpTime, 1);
    
    tmpExpConcs(isnan(tmpExpConcs))=0;
    tmpSDs(tmpExpConcs==0)=0;
    
    vecExpConcs(mm+(1:nExpTime)) = tmpExpConcs;
    vecSDs(mm+(1:nExpTime)) = tmpSDs;
    vecIdMet(mm+(1:nExpTime)) = idMet;
    vecIdTime(mm+(1:nExpTime)) = (1:nExpTime)';
    vecIsComp(mm+(1:nExpTime)) = isComp;
    
    mm = mm + nExpTime;
end
nExpConcs = mm;

vecExpConcStruct.vecExpConcs = vecExpConcs;
vecExpConcStruct.vecSDs = vecSDs;
vecExpConcStruct.vecIdMet = vecIdMet;
vecExpConcStruct.vecIdTime = vecIdTime;
vecExpConcStruct.vecIsComp = vecIsComp;
vecExpConcStruct.nExpConcs = nExpConcs;

end
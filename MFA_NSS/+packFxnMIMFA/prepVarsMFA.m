function varSet = prepVarsMFA170803(model, expData, optionsMFA)

tmpModel = recoverRxnMetInfo(model);

%% ÇÊÇ≠égÇ§ïœêî
varSet.nMets = length(tmpModel.mets);
varSet.nRxns = length(tmpModel.rxns);
varSet.idIndFluxes = tmpModel.idIndFluxes;
varSet.nIndFluxes = length(tmpModel.idIndFluxes);
varSet.idNonIndFluxes = setdiff(1:varSet.nRxns, tmpModel.idIndFluxes);
varSet.nNonIndFluxes = length(varSet.idNonIndFluxes);
varSet.idNetRxns = tmpModel.idNetRxns;
varSet.nNetRxns = length(tmpModel.idNetRxns);
varSet.idNonPoolMets =find(~tmpModel.isPoolMets);
varSet.nNonPoolMets = nnz(~tmpModel.isPoolMets);
varSet.nKnots = optionsMFA.nKnots;
varSet.nCompt = length(expData.comptName);
varSet.idComptCorrParam = optionsMFA.idComptCorrParam;
varSet.nComptCorrParam = length(optionsMFA.idComptCorrParam);
varSet.isEvalConcMets = tmpModel.isEvalConc;
varSet.idEvalConcMets = find(tmpModel.isEvalConc);
varSet.idNonEvalConcMets = find(~tmpModel.isEvalConc(~tmpModel.isPoolMets));
varSet.nEvalConcMets = nnz(tmpModel.isEvalConc);
varSet.nNonEvalConcMets = nnz(~tmpModel.isEvalConc(~tmpModel.isPoolMets));
varSet.isEvalMDVMets = tmpModel.massType == 0;
varSet.idEvalMDVMets = find(tmpModel.massType == 0);
varSet.nEvalMDVMets = nnz(tmpModel.massType == 0);
varSet.idSameEval= tmpModel.idSameEval;


% mergeÇµÇƒï]âøÇ∑ÇÈë„é”ï® (AcCoAÇ∆Ç©)Ç™Ç†ÇÈèÍçáÇÕÅAinnerOptimMetsÇ…ÇÕ1î‘ñ⁄ÇÃÇ›Çä‹ÇﬂÇÈ
idSameEval = varSet.idSameEval;
unqIdSameEval = unique(idSameEval(idSameEval>=1));
TF = idSameEval == 0;
for i = 1 : length(unqIdSameEval)
    id = find(idSameEval==unqIdSameEval(i));
    TF(id(1)) = true;
end

if optionsMFA.isIncludeAllInitConcsInParamMH
    varSet.isInnerOptimMetsInitConcs = false(varSet.nMets,1);
else
    varSet.isInnerOptimMetsInitConcs = ismember(tmpModel.compt, optionsMFA.idComptInnerOptim) & tmpModel.isEvalConc;
end
varSet.isInnerOptimMetsInitConcs = varSet.isInnerOptimMetsInitConcs & TF; % mergeë„é”ï®ÇÃ2î‘ñ⁄à»ç~ÇèúäO

varSet.isInnerOptimMetsConcRates = ismember(tmpModel.compt, optionsMFA.idComptInnerOptim) & tmpModel.isEvalConc;
varSet.isInnerOptimMetsConcRates = varSet.isInnerOptimMetsConcRates & TF; % mergeë„é”ï®ÇÃ2î‘ñ⁄à»ç~ÇèúäO
varSet.idInnerOptimMetsInitConcs = find(varSet.isInnerOptimMetsInitConcs);
varSet.idInnerOptimMetsConcRates =find(varSet.isInnerOptimMetsConcRates);
varSet.nInnerOptimMetsInitConcs = nnz(varSet.isInnerOptimMetsInitConcs);
varSet.nInnerOptimMetsConcRates = nnz(varSet.isInnerOptimMetsConcRates);

varSet.isOuterOptimMetsInitConcs = ~varSet.isInnerOptimMetsInitConcs & ~tmpModel.isPoolMets;
varSet.isOuterOptimMetsConcRates = ~varSet.isInnerOptimMetsConcRates & ~tmpModel.isPoolMets;
varSet.idOuterOptimMetsInitConcs = find(varSet.isOuterOptimMetsInitConcs);
varSet.idOuterOptimMetsConcRates =find(varSet.isOuterOptimMetsConcRates);
varSet.nOuterOptimMetsInitConcs = nnz(varSet.isOuterOptimMetsInitConcs);
varSet.nOuterOptimMetsConcRates = nnz(varSet.isOuterOptimMetsConcRates);


varSet.nExpTime = length(expData.time);

varSet = rmfield(varSet, {'isInnerOptimMetsInitConcs', 'isInnerOptimMetsConcRates','isOuterOptimMetsInitConcs', 'isOuterOptimMetsConcRates'});

end

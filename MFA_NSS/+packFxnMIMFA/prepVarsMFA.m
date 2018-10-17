%% prepare frequently-used variables 
function varSet = prepVarsMFA170803(model, expData, optionsMFA)

tmpModel = recoverRxnMetInfo(model);

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
varSet.nSwitchTimes = optionsMFA.nSwitchTimes;
varSet.isEvalConcMets = tmpModel.isEvalConc;
varSet.idEvalConcMets = find(tmpModel.isEvalConc);
varSet.idNonEvalConcMets = find(~tmpModel.isEvalConc(~tmpModel.isPoolMets));
varSet.nEvalConcMets = nnz(tmpModel.isEvalConc);
varSet.nNonEvalConcMets = nnz(~tmpModel.isEvalConc(~tmpModel.isPoolMets));
varSet.isEvalMDVMets = tmpModel.massType == 0;
varSet.idEvalMDVMets = find(tmpModel.massType == 0);
varSet.nEvalMDVMets = nnz(tmpModel.massType == 0);
varSet.idSameEval= tmpModel.idSameEval;

if optionsMFA.isIncludeAllInitConcsInParamMH
    varSet.isInnerOptimMetsInitConcs = false(varSet.nMets,1);
else
    varSet.isInnerOptimMetsInitConcs =  tmpModel.isEvalConc;
end
varSet.isInnerOptimMetsInitConcs = varSet.isInnerOptimMetsInitConcs;

varSet.isInnerOptimMetsConcRates = tmpModel.isEvalConc;
varSet.isInnerOptimMetsConcRates = varSet.isInnerOptimMetsConcRates ; 
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

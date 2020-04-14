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

if optionsMFA.isUseConcAsParam
    if optionsMFA.isIncludeInitEndConcsInParamMH
        varSet.isQPMetsConcs = false(varSet.nMets,1);
    else
        varSet.isQPMetsConcs = tmpModel.isEvalConc;
    end
    varSet.isQPMetsConcRates = tmpModel.isEvalConc;
else
    if optionsMFA.isIncludeInitEndConcsInParamMH
        varSet.isQPMetsConcs = false(varSet.nMets,1);
    else
        varSet.isQPMetsConcs =  tmpModel.isEvalConc;
    end
    varSet.isQPMetsConcRates = tmpModel.isEvalConc;
end
varSet.idQPMetsConcs = find(varSet.isQPMetsConcs);
varSet.idQPMetsConcRates =find(varSet.isQPMetsConcRates);
varSet.nQPMetsConcs = nnz(varSet.isQPMetsConcs);
varSet.nQPMetsConcRates = nnz(varSet.isQPMetsConcRates);

varSet.isMHMetsConcs = ~varSet.isQPMetsConcs & ~tmpModel.isPoolMets;
varSet.isMHMetsConcRates = ~varSet.isQPMetsConcRates & ~tmpModel.isPoolMets;
varSet.idMHMetsConcs = find(varSet.isMHMetsConcs);
varSet.idMHMetsConcRates =find(varSet.isMHMetsConcRates);
varSet.nMHMetsConcs = nnz(varSet.isMHMetsConcs);
varSet.nMHMetsConcRates = nnz(varSet.isMHMetsConcRates);


varSet.nExpTime = length(expData.time);

varSet = rmfield(varSet, {'isQPMetsConcs', 'isQPMetsConcRates','isMHMetsConcs', 'isMHMetsConcRates'});

end

function [matVar2NotForOpt] = ...
    makeVar2NotForOptim(model, expData, optionsMFA, fullSwitchTimes, varType)

idIndFluxes = optionsMFA.varSet.idIndFluxes ;
idNonPoolMets = optionsMFA.varSet.idNonPoolMets;
nRxns = optionsMFA.varSet.nRxns;
switch varType
    case {'indFlux'}
        matNSwitchTimes = optionsMFA.matNSwitchTimesRxns(idIndFluxes,:);
    case {'concRate'}
        matNSwitchTimes = optionsMFA.matNSwitchTimesMets(idNonPoolMets,:);
    case {'nonIndflux'}
        loc = true(nRxns, 1);
        loc(idIndFluxes) = false;
        matNSwitchTimes = optionsMFA.matNSwitchTimesRxns(loc,:);
    case {'flux'}
        matNSwitchTimes = optionsMFA.matNSwitchTimesRxns;
end

nVar = numel(matNSwitchTimes);
nVarNotForOpt = nnz(~matNSwitchTimes);
matVar2NotForOpt = zeros(nVarNotForOpt, nVar);

rr = 0;
for r = 1 : size(matNSwitchTimes,1)
    if all(matNSwitchTimes(r,:)==true)
        continue
    end
    idFullSwitchTimes = find(~matNSwitchTimes(r,:));
    for k = 1 : length(idFullSwitchTimes)
        idt0 = idFullSwitchTimes(k)-1;
        while true
            t0 = fullSwitchTimes(idt0);
            if matNSwitchTimes(r,idt0)
                break
            end
            idt0 = idt0-1;
        end
        idt2 = idFullSwitchTimes(k)+1;
        while true
            t2 = fullSwitchTimes(idt2);
            if matNSwitchTimes(r,idt2)
                break
            end
            idt2 = idt2+1;
        end
        idt1 = idFullSwitchTimes(k);
        t1 = fullSwitchTimes(idt1);
        
        tmpLoct1 = false(size(matNSwitchTimes));
        tmpLoct1(r,idt1) = true;
        loct1=reshape(tmpLoct1, numel(tmpLoct1),1);
        
        rr = rr +1;
        matVar2NotForOpt(rr, loct1) = 1;
    end
end
matVar2NotForOpt = sparse(matVar2NotForOpt);
end
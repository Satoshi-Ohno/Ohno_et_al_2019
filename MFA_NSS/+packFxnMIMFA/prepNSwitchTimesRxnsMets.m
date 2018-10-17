function     optionsMFA = prepNSwitchTimesRxnsMets(model, expData, optionsMFA)

field2var(optionsMFA.varSet)
if ~optionsMFA.isRxnDepNSwitchTimes
    optionsMFA.matNSwitchTimesRxns = true(nRxns, nSwitchTimes+2);
    optionsMFA.matNSwitchTimesNetRxns = true(nNetRxns, nSwitchTimes+2);
    optionsMFA.matNSwitchTimesMets = true(nMets, nSwitchTimes+2);
    return
end

matNSwitchTimesRxns = true(nRxns, nSwitchTimes+2);
for k = 1 : optionsMFA.nSwitchTimes
    loc = optionsMFA.nSwitchTimesRxns == optionsMFA.patternNSwitchTimes(k);
    matNSwitchTimesRxns(loc, 2+optionsMFA.patternNSwitchTimes(k):end-1) = false;
end
optionsMFA.matNSwitchTimesRxns = matNSwitchTimesRxns;

matNSwitchTimesNetRxns = true(nNetRxns, nSwitchTimes+2);
for k = 1 : optionsMFA.nSwitchTimes
    loc = optionsMFA.nSwitchTimesNetRxns == optionsMFA.patternNSwitchTimes(k);
    matNSwitchTimesNetRxns(loc, 2+optionsMFA.patternNSwitchTimes(k):end-1) = false;
end
optionsMFA.matNSwitchTimesNetRxns = matNSwitchTimesNetRxns;


matNSwitchTimesMets = true(nMets, nSwitchTimes+2);
for k = 1 : optionsMFA.nSwitchTimes
    loc = optionsMFA.nSwitchTimesMets == optionsMFA.patternNSwitchTimes(k);
    matNSwitchTimesMets(loc, 2+optionsMFA.patternNSwitchTimes(k):end-1) = false;
end
optionsMFA.matNSwitchTimesMets = matNSwitchTimesMets;

end

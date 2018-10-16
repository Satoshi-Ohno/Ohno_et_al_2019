function     optionsMFA = prepNKnotsRxnsMets(model, expData, optionsMFA)

field2var(optionsMFA.varSet)
if ~optionsMFA.isRxnSpecificNKnots
    optionsMFA.matNKnotsRxns = true(nRxns, nKnots+2);
    optionsMFA.matNKnotsNetRxns = true(nNetRxns, nKnots+2);
    optionsMFA.matNKnotsMets = true(nMets, nKnots+2);
    return
%     optionsMFA.nKnotsRxns = zeros(nRxns,1) + nKnots;
%     optionsMFA.nKnotsNetRxns = zeros(nNetRxns,1) + nKnots;
%     optionsMFA.nKnotsMets = zeros(nMets,1) + nKnots;
end

matNKnotsRxns = true(nRxns, nKnots+2);
for k = 1 : optionsMFA.nKnots
    loc = optionsMFA.nKnotsRxns == optionsMFA.patternNKnots(k);
    matNKnotsRxns(loc, 2+optionsMFA.patternNKnots(k):end-1) = false;
end
optionsMFA.matNKnotsRxns = matNKnotsRxns;

matNKnotsNetRxns = true(nNetRxns, nKnots+2);
for k = 1 : optionsMFA.nKnots
    loc = optionsMFA.nKnotsNetRxns == optionsMFA.patternNKnots(k);
    matNKnotsNetRxns(loc, 2+optionsMFA.patternNKnots(k):end-1) = false;
end
optionsMFA.matNKnotsNetRxns = matNKnotsNetRxns;


matNKnotsMets = true(nMets, nKnots+2);
for k = 1 : optionsMFA.nKnots
    loc = optionsMFA.nKnotsMets == optionsMFA.patternNKnots(k);
    matNKnotsMets(loc, 2+optionsMFA.patternNKnots(k):end-1) = false;
end
optionsMFA.matNKnotsMets = matNKnotsMets;

end

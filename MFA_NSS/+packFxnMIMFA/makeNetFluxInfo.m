%% Make net flux information
function model = makeNetFluxInfo170821(model, optionsMFA)

[nMets, nRxns] = size(model.S);
revSets=model.rxnInfo.revSets;

%% Identification of reversible reactions
idRxns = (1:nRxns)';
idNonRevRxns = idRxns(revSets == 0);
[~,idUnqRevRxns] = unique(revSets);
idNetRxns = sort([idNonRevRxns;idUnqRevRxns(2:end)]);
nNetRxns = length(idNetRxns);

%%  Make reaction list
netS = model.S(:, idNetRxns);
isRev = revSets(idNetRxns) > 0;
netRxns = Smatrix2rxns(netS, model.mets, isRev);
netRxnNames = model.rxnNames(idNetRxns);

%% Make matrix to convert raw flux to net flux
corrRatio = zeros(nNetRxns,1)+1;
matRaw2Net = zeros(nNetRxns, nRxns+1);

matRaw2Net(1:nNetRxns, idNetRxns) =...
    matRaw2Net(1:nNetRxns, idNetRxns)+ eye(nNetRxns);
matRaw2Net(1:nNetRxns, idNetRxns+1) = ...
    matRaw2Net(1:nNetRxns, idNetRxns+1)-  eye(nNetRxns);
matRaw2Net(ismember(idNetRxns, idNonRevRxns), idNonRevRxns+1) = ...
    matRaw2Net(ismember(idNetRxns, idNonRevRxns), idNonRevRxns+1) + eye(length(idNonRevRxns));
matRaw2Net = matRaw2Net(:,1:nRxns);
matRaw2Net = matRaw2Net .* repmat(corrRatio, 1, nRxns); % 補正リストに従ってフラックスを修正

%% 
model.netRxns = netRxns;
model.netRxnNames = netRxnNames;
model.idNetRxns = idNetRxns;
model.matRaw2Net = matRaw2Net;

end

%% Make reaction formula from stoichiometric matrix
function rxns = Smatrix2rxns(S, mets, isRev)

[nMets, nRxns] = size(S);
rxns = cell(nRxns,1);
for r = 1 : nRxns
    tmpRxn = [];
    idMet = [find(S(:,r)<0) ; find(S(:,r)>0)];
%     S(idMet,r)
    
    if isRev(r)
        arrowChar = '<=>';
    else
        arrowChar = '-->';
    end
    for m = 1 : length(idMet)
        if m == 1
            tmpRxn = mets{idMet(m)};
        else
            if S(idMet(m-1),r) < 0 && S(idMet(m),r) > 0
                tmpRxn = [tmpRxn ' ' arrowChar ' ' mets{idMet(m)}];
            else
                tmpRxn = [tmpRxn ' + ' mets{idMet(m)}];
            end
        end
        
    end
    rxns{r} = tmpRxn;
end

end




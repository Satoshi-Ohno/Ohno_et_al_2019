%% Make transformation matrix (concConst + concRate -> conc)
function matConcConstRate2Conc = ...
    makeMatConcConstRate2Conc(model, expData, optionsMFA, fullSwitchTimes, timeConstr, idMetsAll, isEarlierTI)

nNonPoolMets = optionsMFA.varSet.nNonPoolMets;
nSwitchTimes = optionsMFA.varSet.nSwitchTimes;

idTimeInterval = zeros(1,length(timeConstr));
switchTimeOrder = 1 : nSwitchTimes+1;
for k = switchTimeOrder
    loc = timeConstr>=fullSwitchTimes(k) & timeConstr<=fullSwitchTimes(k+1);
    idTimeInterval(loc) = k;
end


matConcConstRate2Conc = sparse(...
    length(timeConstr)*length(idMetsAll),...
    nNonPoolMets*((nSwitchTimes+1) + (nSwitchTimes+2))...
    );
kk = 0;
for k = 1 : nSwitchTimes+1
    tmpMatCoef = zeros(3,3);
    tmpMatCoef(1,1) = fullSwitchTimes(k+1)-fullSwitchTimes(k);
    tmpMatCoef(2,2) = fullSwitchTimes(k+1);
    tmpMatCoef(2,3) = -fullSwitchTimes(k);
    tmpMatCoef(3,2) = -1;
    tmpMatCoef(3,3) = 1;
    
    loc = idTimeInterval == k;
    tmpTime = timeConstr(loc);
    if isempty(tmpTime)
        continue
    end
    tmpMatTime = [ones(length(tmpTime), 1),  tmpTime',  1/2*tmpTime'.^2 ];
    
    tmpMatTimeCoef = zeros(size(tmpMatTime,1)*length(idMetsAll), size(tmpMatTime,2)*length(idMetsAll));
    for m = 1 : length(idMetsAll)
        row = size(tmpMatTime,1)*(m-1) + (1:size(tmpMatTime,1));
        col = size(tmpMatTime,2)*(m-1) + (1:size(tmpMatTime,2));
        tmpMatTimeCoef(row, col) = tmpMatTimeCoef(row, col) + tmpMatTime * tmpMatCoef;
    end
    
    
    row = repmat(length(timeConstr)*((1:length(idMetsAll))-1)', 1, length(tmpTime)) ...
        + repmat((1:length(tmpTime)), length(idMetsAll), 1) ...
        + kk;
    row = reshape(row', numel(row), 1);
    
    colConcConst = nNonPoolMets*(k-1)+idMetsAll;
    colConcRate1 = nNonPoolMets*(nSwitchTimes+1)+nNonPoolMets*(k-1)+idMetsAll;
    colConcRate2 = nNonPoolMets*(nSwitchTimes+1)+nNonPoolMets*(k)+idMetsAll;
    col = [colConcConst, colConcRate1, colConcRate2];
    col = reshape(col', numel(col), 1);
    
    matConcConstRate2Conc(row, col) = matConcConstRate2Conc(row, col) + ...
        1/(fullSwitchTimes(k+1)-fullSwitchTimes(k)) * sparse(tmpMatTimeCoef);
    
    kk = kk + length(tmpTime);
end


end
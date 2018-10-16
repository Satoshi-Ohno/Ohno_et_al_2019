%% 天然同位体の逆補正行列の作成
function expData = createCorrMatNatIsotope170105(expData, preExpData, model)

nMets = length(model.mets);

%% 天然同位体補正
% 誤差の天然同位体補正はまだ。
% 
%天然同位体比率(de LAETER, J. R., et al.,2003, van winden (2002), BnB)
% representative
natIsotopeRatio.C = [98.93, 1.07]/100;
% natIsotopeRatio.C = [1, 0]/100;
natIsotopeRatio.N = [99.632, 0.368]/100;
natIsotopeRatio.O = [99.757, 0.038, 0.205]/100;
natIsotopeRatio.H = [99.985, 0.015]/100;
% 比率低め
% natIsotopeRatio.C = [99.037, 0.963]/100;
% natIsotopeRatio.N = [99.654, 0.346]/100;
% natIsotopeRatio.O = [99.775, 0.037, 0.188]/100;
% natIsotopeRatio.H = [99.997, 0.003]/100;

corrMat =  cell(size(model.metInfo.metKeggIds));
% preMDVs = MDVs;
% preMDVsSD = MDVsSD;
for i  = 1 : numel(model.metInfo.metKeggIds)
    if isempty(model.metInfo.metKeggIds{i})
        continue
    end
    if rem(i, nMets) == 0
        if ismember(model.metInfo.massType(nMets), [1, 2])
%         if ismember(model.massType(nMets), [1, 2,-1])
            continue
        end
    else
        if ismember(model.metInfo.massType(rem(i, nMets)), [1, 2])
%         if ismember(model.massType(rem(i, nMets)), [1, 2,-1])
            continue
        end
    end
    nMassType = sum(model.metInfo.Ematrix{i})+1;  %質量数は何パターンありうるか
%     tempMDV = zeros(nMassType,length(expData.time));
%     tempMDV(1:model.nCarbonMets(i)+1,:) = preMDVs{i};
    preCorrMat = eye(nMassType);
    for j = 1 : 4
        tmpCorrMat = zeros(nMassType);
        switch j
            case 1 % C
                tmpCorrMat(1:end,1:end-0) = tmpCorrMat(1:end,1:end-0) + eye(nMassType-0) * natIsotopeRatio.C(1);
                tmpCorrMat(2:end,1:end-1) = tmpCorrMat(2:end,1:end-1) + eye(nMassType-1) * natIsotopeRatio.C(2);
            case 2 % N
                tmpCorrMat(1:end,1:end-0) = tmpCorrMat(1:end,1:end-0) + eye(nMassType-0) * natIsotopeRatio.N(1);
                tmpCorrMat(2:end,1:end-1) = tmpCorrMat(2:end,1:end-1) + eye(nMassType-1) * natIsotopeRatio.N(2);
            case 3 % O
                tmpCorrMat(1:end,1:end-0) = tmpCorrMat(1:end,1:end-0) + eye(nMassType-0) * natIsotopeRatio.O(1);
                tmpCorrMat(2:end,1:end-1) = tmpCorrMat(2:end,1:end-1) + eye(nMassType-1) * natIsotopeRatio.O(2);
                tmpCorrMat(3:end,1:end-2) = tmpCorrMat(3:end,1:end-2) + eye(nMassType-2) * natIsotopeRatio.O(3);
            case 4 % H
                tmpCorrMat(1:end,1:end-0) =  tmpCorrMat(1:end,1:end-0) + eye(nMassType-0) * natIsotopeRatio.H(1);
                tmpCorrMat(2:end,1:end-1) =  tmpCorrMat(2:end,1:end-1) + eye(nMassType-1) * natIsotopeRatio.H(2);
        end
        nAtom = model.metInfo.Ematrix{i}(j);
        if j == 1 %炭素のときは骨格炭素をのぞく
            idMet = rem(i, nMets);
            if idMet==0
                idMet = nMets;
            end
            nAtom = nAtom - expData.nCarbonMets(idMet);
        end
        tmpCorrMat = tmpCorrMat^nAtom;
        preCorrMat =  preCorrMat * tmpCorrMat;
%         corrMat{i, j} = sparse(tmpCorrMat);
%         tempMDV = tmpCorrMat\tempMDV;
    end
    corrMat{i} = preCorrMat;
%     corrMat{i} = inv(preCorrMat);
end

expData.nat13CRatio = natIsotopeRatio.C;
expData.corrMat = corrMat;



end

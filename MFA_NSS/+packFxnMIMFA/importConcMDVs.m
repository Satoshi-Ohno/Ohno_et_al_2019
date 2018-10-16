%% 代謝物濃度をインポートを作成
function expData = importConcMDVs160413(expData, preExpData, model, optionsMFA)
%importConcs 実験データより13C標識割合をあらわすMDV (Mass Distribution Vector)を作成
%
% model = importConcs(model, expData, preExpData)
%
% INPUT
%   expData      作成中の実験データ
%       mets                        代謝物名
%       nCarbonMets           代謝物の炭素数
%   preExpData  実験データ
%       cpdName     代謝物名
%       timeCourse
%           avg      濃度
%           sem     標準誤差
%
% OUTPUT
%   expData     実験データ
%       concMDVs          代謝物濃度
%       concMDVsSD      代謝物濃度の標準偏差
%
mets = expData.mets;
nMets = length(mets);
time = expData.time;
nExpTime = length(time);
nRepThreshold = optionsMFA.nRepThreshold ;
outCalib = model.outCalib;

%% conc×MDVのタイムコース

concMDVs = cell(nMets,1);
concMDVsSD = cell(nMets,1);
concMDVsAll = cell(nMets,1);
% minConcMDVsSD = 1;
for m = 1 : nMets
    nCarbon = expData.nCarbonMets(m);
    nMDV = nCarbon+1;
    tmpConcMDVs = zeros(nMDV, nExpTime);
    tmpConcMDVsSD = zeros(nMDV, nExpTime);
    tmpConcMDVsAll = zeros(nMDV, nExpTime,5);
    tmpCompt = expData.comptName(model.metInfo.compt(m));
    switch model.metInfo.massType(m)
        case 1
            tmpConcMDVs(end,:) = 10^6;
        case 2
            tmpConcMDVs(1,:) = 10^6;
        case {0, 3, 4, 5}
            %             cellTmpConcMDVsAll = cell(1,3);
            %             tmpConcMDVsSD = [];
            
            nRep = 0;
            for j = 1 : size(expData.expMets,2)
                if isempty(expData.expMets{m,j}) % 180830追加
                    break
                end
                expMetId = strcmp(expData.expMets(m,j), {preExpData.cpdName}) ...
                    & strcmp(tmpCompt, {preExpData.compartment});
                nRep = max([nRep, size(preExpData(expMetId).timeCourse.concMDVsAll,3)]);
            end
            tmpConcMDVsAll = zeros(nMDV, nExpTime, nRep);
            for j = 1 : size(expData.expMets,2)
                if isempty(expData.expMets{m,j})
                    break
                end
                expMetId = strcmp(expData.expMets(m,j), {preExpData.cpdName}) ...
                    & strcmp(tmpCompt, {preExpData.compartment});
                %                 expMetId = strcmp(expData.expMets(i,j), {preExpData.cpdName}) ...
                %                     & model.isExtracellular(i)==[preExpData.isExtracellular];
                tmp = preExpData(expMetId).timeCourse.concMDVsAll(1:nMDV,:,:);
                tmp(isnan(tmp))=0;
                tmpConcMDVsAll = tmpConcMDVsAll + tmp;
            end
            tmpConcMDVsAll(tmpConcMDVsAll==0) = nan;
            %             tmpConcMDVsAll = zeros(nMDV,nTime,3);
            %             for n = 1 : 3
            %                 tmpConcMDVsAll(:,:,n) = sum(cellTmpConcMDVsAll{n},3);
            %             end
            
            % out of calibrationへの対応
            idOutCalib= find(m==model.outCalib.idMets);
            if ~isempty(idOutCalib)
                for j = 1 : length(idOutCalib)
                    idExp = outCalib.expIds(idOutCalib(j));
                    idTime = find(outCalib.time(idOutCalib(j))==time);
                    tmpConcMDVsAll(:, idTime, idExp) = nan;
                end
            end
            
            % 平均値とSDの算出
            tmpConcMDVs = nan(nMDV, nExpTime);
            tmpConcMDVsSD = nan(nMDV, nExpTime);
            for i = 1 : nMDV
                for t = 1 : nExpTime
                    locIsNan = squeeze(isnan(tmpConcMDVsAll(i,t, :)));
                    if nnz(~locIsNan) >= nRepThreshold
                        %             if nnz(concMDVsAll(i,t,:)>0)==3
                        tmpConcMDVs(i,t) = mean(tmpConcMDVsAll(i,t,~locIsNan));
                        if optionsMFA.isUseInputCV
                            tmpConcMDVsSD(i,t) = tmpConcMDVs(i,t)*optionsMFA.inputConcsCV;
                        else
                            tmpConcMDVsSD(i,t) = std(tmpConcMDVsAll(i,t,~locIsNan));
                            tmpConcMDVsSD(i,t) = max([...
                                tmpConcMDVsSD(i,t), tmpConcMDVs(i,t)*optionsMFA.minConcsCV]); % 最小のSDを設定
                        end
                    end
                end
            end
    end
    concMDVs{m} = tmpConcMDVs;
    concMDVsSD{m} = tmpConcMDVsSD;
    concMDVsAll{m} = tmpConcMDVsAll;
    
end

% % 実験データがない代謝物について、その代謝物の最小の誤差を与える。
% for i = 1 : nMets
%     concMDVsSD{i}(concMDVsSD{i} == 0) = minConcMDVsSD;
% end

expData.concMDVs = concMDVs;
expData.concMDVsSD = concMDVsSD;
expData.concMDVsAll = concMDVsAll;
end

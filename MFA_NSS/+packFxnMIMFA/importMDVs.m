%% MDVをインポート
function expData = importMDVs160817(expData, preExpData, model, optionsMFA)
%importMDVs 実験データより13C標識割合をあらわすMDV (Mass Distribution Vector)を作成
% 
% model = importMDVs(model, expData, preExpData)
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
%       expMets     実験データの代謝物名
%       MDVs        13C標識割合
%       MDVsSD      13C標識割合の標準偏差
% 
% 2016.4.13  MDV算出時に、nMDVまでの総和が1になるように補正

nMets = length(model.mets);
time = expData.time;
nExpTime = length(time);
outCalib = model.outCalib;
nRepThreshold = optionsMFA.nRepThreshold ;
minMDVsSD = optionsMFA.minMDVsSD;

%% MDVのタイムコース

MDVs = cell(nMets,1);
MDVsSD = cell(nMets,1);
MDVsAll = cell(nMets,1);
% minMDVsSD = 1;
for m = 1 : nMets
    nMDV = expData.nCarbonMets(m)+1;
    tmpMDVs = zeros(nMDV, nExpTime);
    tmpMDVsSD = zeros(nMDV, nExpTime);
    nRep = size(expData.concMDVsAll{m},3);
    tmpMDVsAll = zeros(nMDV, nExpTime,nRep);
%     tmpCompartment = expData.compartmentName(model.compartment(i));
    switch model.metInfo.massType(m)
        case 1
            tmpMDVs(end,:) = 1;
        case 2
            tmpMDVs(1,:) = 1;
        case {0,3, 4, 5}
%             tmpMDVsAll = zeros(nMDV, nExpTime, nRep);
            for n = 1 : nRep
                tmpConcMDVs = expData.concMDVsAll{m}(:,:,n);
                tmpConcMDVs(isnan(tmpConcMDVs))=0;
                tmpConcs = sum(tmpConcMDVs(1:nMDV,:),1);
                tmpMDVsAll(:,:,n) = tmpConcMDVs ./ repmat(tmpConcs, nMDV,1);
%                 tmpMDVsN3(:,:,n) = expData.concMDVsN3{i}(:,:,n) ./ repmat(expData.concsN3{i}(n,:), nMDV,1);
            end
            tmpMDVsAll(tmpMDVsAll==0) = nan;
            
            % out of calibrationへの対応
            idOutCalib= find(m==model.outCalib.idMets);
            if ~isempty(idOutCalib)
                for j = 1 : length(idOutCalib)
                    idExp = outCalib.expIds(idOutCalib(j));
                    idTime = find(outCalib.time(idOutCalib(j))==time);
                    tmpMDVsAll(:, idTime, idExp) = nan;
                end
            end

            % 平均値とSDの算出
            tmpMDVs = nan(nMDV, nExpTime);
            tmpMDVsSD = nan(nMDV, nExpTime);
            for i = 1 : nMDV
                for t = 1 : nExpTime
                    locIsNan = squeeze(isnan(tmpMDVsAll(i,t, :)));
                    if nnz(~locIsNan) >= nRepThreshold
                        tmpMDVs(i,t) = mean(tmpMDVsAll(i,t,~locIsNan));
                        if optionsMFA.isUseInputCV
                            tmpMDVsSD(i,t) = tmpMDVs(i,t)*optionsMFA.inputMDVsCV;
                        else
                            tmpMDVsSD(i,t) = std(tmpMDVsAll(i,t,~locIsNan));
                            tmpMDVsSD(i,t) = max([tmpMDVsSD(i,t), optionsMFA.minMDVsSD]); % 最小のSDを設定
                        end
                    end
                end
            end
    end
    tmpMDVs(tmpMDVsSD==0) = nan;
    tmpMDVsSD(tmpMDVsSD==0) = nan;
    MDVs{m} = tmpMDVs;
    MDVsSD{m} = tmpMDVsSD;
    MDVsAll{m} = tmpMDVsAll;
    
end

expData.MDVs = MDVs;
expData.MDVsSD = MDVsSD;
expData.MDVsAll = MDVsAll;

end





%% 代謝物濃度をインポートを作成
function expData = importConcs160413(expData, preExpData, model, optionsMFA)
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
%       Concs          代謝物濃度
%       ConcsSD      代謝物濃度の標準偏差
% 

mets = expData.mets;
nMets = length(mets);
time = expData.time;
nExpTime = length(time);
outCalib = model.outCalib;
nRepThreshold = optionsMFA.nRepThreshold ;

% 代謝物濃度のタイムコース
concsAll = cell(size(model.metInfo.metKeggIds));
concsSDOri = cell(size(model.metInfo.metKeggIds));
concs = zeros(nMets,length(time));
concsSD = zeros(nMets,length(time));
concsAll = cell(nMets, 1);
for i = 1 : nMets
    tmpConcs = zeros(1, length(time));
    tmpConcsSD = zeros(1, length(time));
    tmpConcsAll = zeros(1, length(time),5);
    tmpCompt = expData.comptName(model.metInfo.compt(i));
    switch model.metInfo.isEvalConc(i)
%     switch model.massType(i)
        case {0}
%         case {1,2}
            tmpConcs(1,:) = 10^6;
        case {1}
%         case {0,3,4}
%             timeCourse = zeros(1,nTime, nRep);
            nRep=0;
            for j = 1 : size(expData.expMets,2)
                if isempty(expData.expMets{i,j}) % 180830追加
                    break
                end
                expMetId = strcmp(expData.expMets(i,j), {preExpData.cpdName}) ...
                    & strcmp(tmpCompt, {preExpData.compartment});
                nRep = max([nRep, size(preExpData(expMetId).timeCourse.concsAll,3)]);
            end
            tmpConcsAll = zeros(1,nExpTime, nRep);
            for j = 1 : size(expData.expMets,2)
                if isempty(expData.expMets{i,j})
                    break
                end
                expMetId = strcmp(expData.expMets(i,j), {preExpData.cpdName}) ...
                    & strcmp(tmpCompt, {preExpData.compartment});
%                 expMetId = strcmp(expData.expMets(i,j), {preExpData.cpdName}) ...
%                     & model.isExtracellular(i)==[preExpData.isExtracellular];
                if nnz(expMetId) == 0
                    error('Metabolite name is not in the exp data file.')
                end
                if nnz(expMetId) ~= 1
                    error('Metabolite name is not unique.')
                end
                tmp = preExpData(expMetId).timeCourse.concsAll;
                tmp(isnan(tmp))=0;
                tmpConcsAll= tmpConcsAll+tmp;
%                 timeCourse = cat(3,timeCourse, preExpData(expMetId).timeCourse.concsAll);
            end
            tmpConcsAll(tmpConcsAll==0) = nan;
%             tmpConcsAll = timeCourse;
%             tmpConcsAll = sum(timeCourse, 3);
            
            %out of calibrationへの対応
            idOutCalib= find(i==model.outCalib.idMets);
             if ~isempty(idOutCalib)
                 for j = 1 : length(idOutCalib)
                     idExp = outCalib.expIds(idOutCalib(j));
                     idTime = find(outCalib.time(idOutCalib(j))==time);
                     tmpConcsAll(idExp, idTime) = nan;
                 end
             end
             
             % 平均値とSDの算出
             tmpConcs = nan(1,nExpTime);
             tmpConcsSD = nan(1,nExpTime);
             for t = 1 : nExpTime
                 locIsNan = squeeze(isnan(tmpConcsAll(1,t, :)));
                 if nnz(~locIsNan) >= nRepThreshold
                     tmpConcs(t) = mean(tmpConcsAll(1,t,~locIsNan));
                     if optionsMFA.isUseInputCV
                         tmpConcsSD(t) = tmpConcs(t)*optionsMFA.inputConcsCV;
                     else
                         tmpConcsSD(t) = std(tmpConcsAll(1,t,~locIsNan));
                         tmpConcsSD(t) = max([...
                             tmpConcsSD(t), tmpConcs(t)*optionsMFA.minConcsCV]); % 最小のSDを設定
                     end
                 end
             end
%              tmpConcs = nan(1, nTime);
%             tmpConcsSD = nan(1, nTime);
%             for t = 1 : nTime
%                 if nnz(tmpConcsAll(1,t,:)>0)>=3
%                     isNanExp = isnan(tmpConcsAll(1,t,:)) | tmpConcsAll(1,t,:)==0;
%                     tmpConcs(t) = mean(tmpConcsAll(1,t,~isNanExp));
%                     if optionsMFA.isUseInputCV
%                         tmpConcsSD(t) = tmpConcs(t)*optionsMFA.inputConcsCV;
%                     else
%                         tmpConcsSD(t) = std(tmpConcsAll(1,t,~isNanExp));
%                     end
%                 end
%             end
    end
    concs(i,:) = tmpConcs;
    concsSD(i,:) = tmpConcsSD;
    concsAll{i} = tmpConcsAll;
end
expData.concs = concs;
expData.concsSD = concsSD;
expData.concsAll = concsAll;

end



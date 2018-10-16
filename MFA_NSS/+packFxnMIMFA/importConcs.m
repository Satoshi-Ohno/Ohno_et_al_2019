%% ��ӕ��Z�x���C���|�[�g���쐬
function expData = importConcs160413(expData, preExpData, model, optionsMFA)
%importConcs �����f�[�^���13C�W������������킷MDV (Mass Distribution Vector)���쐬
% 
% model = importConcs(model, expData, preExpData)
% 
% INPUT
%   expData      �쐬���̎����f�[�^
%       mets                        ��ӕ���
%       nCarbonMets           ��ӕ��̒Y�f��
%   preExpData  �����f�[�^
%       cpdName     ��ӕ���
%       timeCourse  
%           avg      �Z�x
%           sem     �W���덷
% 
% OUTPUT
%   expData     �����f�[�^
%       Concs          ��ӕ��Z�x
%       ConcsSD      ��ӕ��Z�x�̕W���΍�
% 

mets = expData.mets;
nMets = length(mets);
time = expData.time;
nExpTime = length(time);
outCalib = model.outCalib;
nRepThreshold = optionsMFA.nRepThreshold ;

% ��ӕ��Z�x�̃^�C���R�[�X
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
                if isempty(expData.expMets{i,j}) % 180830�ǉ�
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
            
            %out of calibration�ւ̑Ή�
            idOutCalib= find(i==model.outCalib.idMets);
             if ~isempty(idOutCalib)
                 for j = 1 : length(idOutCalib)
                     idExp = outCalib.expIds(idOutCalib(j));
                     idTime = find(outCalib.time(idOutCalib(j))==time);
                     tmpConcsAll(idExp, idTime) = nan;
                 end
             end
             
             % ���ϒl��SD�̎Z�o
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
                             tmpConcsSD(t), tmpConcs(t)*optionsMFA.minConcsCV]); % �ŏ���SD��ݒ�
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



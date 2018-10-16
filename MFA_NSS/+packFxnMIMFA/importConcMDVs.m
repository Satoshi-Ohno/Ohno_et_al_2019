%% ��ӕ��Z�x���C���|�[�g���쐬
function expData = importConcMDVs160413(expData, preExpData, model, optionsMFA)
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
%       concMDVs          ��ӕ��Z�x
%       concMDVsSD      ��ӕ��Z�x�̕W���΍�
%
mets = expData.mets;
nMets = length(mets);
time = expData.time;
nExpTime = length(time);
nRepThreshold = optionsMFA.nRepThreshold ;
outCalib = model.outCalib;

%% conc�~MDV�̃^�C���R�[�X

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
                if isempty(expData.expMets{m,j}) % 180830�ǉ�
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
            
            % out of calibration�ւ̑Ή�
            idOutCalib= find(m==model.outCalib.idMets);
            if ~isempty(idOutCalib)
                for j = 1 : length(idOutCalib)
                    idExp = outCalib.expIds(idOutCalib(j));
                    idTime = find(outCalib.time(idOutCalib(j))==time);
                    tmpConcMDVsAll(:, idTime, idExp) = nan;
                end
            end
            
            % ���ϒl��SD�̎Z�o
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
                                tmpConcMDVsSD(i,t), tmpConcMDVs(i,t)*optionsMFA.minConcsCV]); % �ŏ���SD��ݒ�
                        end
                    end
                end
            end
    end
    concMDVs{m} = tmpConcMDVs;
    concMDVsSD{m} = tmpConcMDVsSD;
    concMDVsAll{m} = tmpConcMDVsAll;
    
end

% % �����f�[�^���Ȃ���ӕ��ɂ��āA���̑�ӕ��̍ŏ��̌덷��^����B
% for i = 1 : nMets
%     concMDVsSD{i}(concMDVsSD{i} == 0) = minConcMDVsSD;
% end

expData.concMDVs = concMDVs;
expData.concMDVsSD = concMDVsSD;
expData.concMDVsAll = concMDVsAll;
end

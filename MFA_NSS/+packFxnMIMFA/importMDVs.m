%% MDV���C���|�[�g
function expData = importMDVs160817(expData, preExpData, model, optionsMFA)
%importMDVs �����f�[�^���13C�W������������킷MDV (Mass Distribution Vector)���쐬
% 
% model = importMDVs(model, expData, preExpData)
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
%       expMets     �����f�[�^�̑�ӕ���
%       MDVs        13C�W������
%       MDVsSD      13C�W�������̕W���΍�
% 
% 2016.4.13  MDV�Z�o���ɁAnMDV�܂ł̑��a��1�ɂȂ�悤�ɕ␳

nMets = length(model.mets);
time = expData.time;
nExpTime = length(time);
outCalib = model.outCalib;
nRepThreshold = optionsMFA.nRepThreshold ;
minMDVsSD = optionsMFA.minMDVsSD;

%% MDV�̃^�C���R�[�X

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
            
            % out of calibration�ւ̑Ή�
            idOutCalib= find(m==model.outCalib.idMets);
            if ~isempty(idOutCalib)
                for j = 1 : length(idOutCalib)
                    idExp = outCalib.expIds(idOutCalib(j));
                    idTime = find(outCalib.time(idOutCalib(j))==time);
                    tmpMDVsAll(:, idTime, idExp) = nan;
                end
            end

            % ���ϒl��SD�̎Z�o
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
                            tmpMDVsSD(i,t) = max([tmpMDVsSD(i,t), optionsMFA.minMDVsSD]); % �ŏ���SD��ݒ�
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





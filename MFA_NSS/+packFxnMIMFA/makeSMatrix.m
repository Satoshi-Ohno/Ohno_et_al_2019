%% ���w�ʘ_�s��쐬
function model = makeSMatrix170301(model, optionsMFA)
%makeSMatrix    ���w�ʘ_�s��쐬
% 
% model = makeSMatrix(model)
% 
% INPUT
%   model       ��Ӄ��f��
%       rxns        ��Ӕ��������X�g
%       carbonTrans ���i�Y�f�̈ړ����X�g
% 
% OUTPUT
%   model       ��Ӄ��f��
%       mets            ��ӕ����X�g
%       S               ���w�ʘ_�s��
%       arrowChar       ��������킷������
%       nCarbonMets     ��ӕ��̒Y�f�����X�g
% 

tmpModel = recoverRxnMetInfo(model);

%% ��ӕ����X�g�̍쐬

% ���������X�y�[�X�ŋ�؂�
arrowChar = '-->';
rxns = tmpModel.rxns;
cellRxns = cell(length(rxns),20);
for j = 1 : length(rxns)
    splitRxn = strsplit(rxns{j},' ');
    cellRxns(j,1:length(splitRxn)) = splitRxn;
end

% ��ӕ����X�g���쐬
isEmptyCell = cellfun('isempty',cellRxns);
mets = unique(cellRxns(~isEmptyCell));
mets = mets(~strcmp(mets, {arrowChar}));
mets = mets(~strcmp(mets, {'+'}));
mets = mets(~strncmp(mets, {'('}, 1));
% if length(mets) ~= length(model.mets)
if ~isempty(setdiff(mets, tmpModel.mets)) || ~isempty(setdiff(tmpModel.mets, mets)) 
    setdiff(mets, tmpModel.mets)
    setdiff(tmpModel.mets, mets)
    error('Metabolite list and/or mass balance equation is not correct.')
%     mets = model.mets;
end

%% ���i�Y�f�̈ړ����X�g�̐���

nRxns = length(rxns);
nMets = length(tmpModel.mets);
carbonTrans = tmpModel.rxnInfo.carbonTransitions;

cellCarbonTrans = cell(nRxns,20);
arrowIds = zeros(nRxns,1);
for j = 1 : nRxns
    splitCarbonTrans = strsplit(carbonTrans{j},' ');
    cellCarbonTrans(j,1:length(splitCarbonTrans)) = splitCarbonTrans;
    arrowIds(j) = find(strcmp(splitCarbonTrans,{arrowChar}));
end


%% ���w�ʘ_�s��ƒY�f�����X�g���쐬
% S: ���w�ʘ_�s��
% Sorder: ����j�ɂ����đ�ӕ�i�����Ԗڂ̊or��������
% Ssubs: ���w�ʘ_�s��̂����A��̌W���݂̂����������́B
% Sprod: ���w�ʘ_�s��̂����A�������̌W���݂̂����������́B

S = zeros(nMets, nRxns);
Sorder = zeros(nMets, nRxns);
Ssubs = zeros(nMets, nRxns);
Sprod = zeros(nMets, nRxns);
idSynsRxns = cell(nMets, 1);
% nCarbonMets = zeros(nMets,1);
for j = 1 : length(rxns)
    tmpOrder = 1;
    posNeg = -1;
    for i = 1 : size(cellRxns,2);
        coef = 1;
        if isempty(cellRxns{j,i})
            break
        elseif strcmp(cellRxns(j,i), {arrowChar})
            posNeg = 1;
            tmpOrder = 1;
            continue
        end
        isMet = strcmp(cellRxns(j,i),tmpModel.mets);
        if any(isMet)
            if i >=2 && strncmp(cellRxns(j,i-1), {'('},1)
                coef = str2double(cellRxns{j,i-1}(2:end-1));
            end
            S(isMet, j) = S(isMet, j) + posNeg*coef;
            Sorder(isMet, j) = Sorder(isMet, j) + posNeg*coef*tmpOrder;
            if posNeg < 0
            Ssubs(isMet, j) = Ssubs(isMet, j) + posNeg*coef;
            else
            Sprod(isMet, j) = Sprod(isMet, j) + posNeg*coef;
            end
            tmpOrder = tmpOrder+1;
        end
    end
end

model.S = sparse(S);
model.Sorder = sparse(Sorder);
model.Ssubs = sparse(Ssubs);
model.Sprod = sparse(Sprod);
model.arrowChar = arrowChar;
end



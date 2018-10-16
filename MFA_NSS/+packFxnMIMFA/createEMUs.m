%% ICM���쐬
function model = createEMUs161129(model, optionsMFA)
%createEMUs EMU�l�b�g���[�N���\�z (Isotopomer Converted Matrix)���쐬
% 
% model = createEMUs(model)
% 
% INPUT
%   model       ��Ӄ��f��
%       carbonTransitions   ���i�Y�f�̈ړ���\����
%       S                   ���w�ʘ_�s��
%       nCarbonMets         ��ӕ��̒Y�f�����X�g
% 
% OUTPUT
%   model       ��Ӄ��f��
%       emuNetwork        EMU�̃l�b�g���[�N
% 

tmpModel = recoverRxnMetInfo(model);


[~, nRxns] = size(tmpModel.S);

%% ���������X�y�[�X�ŋ�؂�

rxns = tmpModel.rxns;
carbonTrans = tmpModel.carbonTransitions;

cellRxns = cell(length(rxns),20);
cellCTrans = cell(nRxns,20);
arrowIds = zeros(nRxns,1);
for r = 1 : length(rxns)
    splitRxn = strsplit(rxns{r},' ');
    cellRxns(r,1:length(splitRxn)) = splitRxn;
    
    splitCarbonTrans = strsplit(carbonTrans{r},' ');
    cellCTrans(r,1:length(splitCarbonTrans)) = splitCarbonTrans;
    
    arrowIds(r) = find(strcmp(splitRxn,{tmpModel.arrowChar}));
end

% ���L�t�B�[���h�͉��ō폜
tmpModel.cellRxns = cellRxns;
tmpModel.cellCTrans = cellCTrans;
tmpModel.arrowIds = arrowIds;

%% ���͂����ӕ��ɂ���EMU���`

idNonPoolMets = find(tmpModel.massType<=0);
% idMesMet = find(model.massType==0);
for m = 1 : length(idNonPoolMets)
    mesEmuList(m).name = [];
    mesEmuList(m).met = tmpModel.mets{idNonPoolMets(m)};
    mesEmuList(m).isCInEmu = true(1,tmpModel.nCarbonMets(idNonPoolMets(m)));
    mesEmuList(m).name = char(defineEmu(mesEmuList(m).met, mesEmuList(m).isCInEmu));
    mesEmuList(m).size = tmpModel.nCarbonMets(idNonPoolMets(m));
    mesEmuList(m).isSmallerSize = false;
    mesEmuList(m).isKnownEmu = false;
end
mesEmuSizes = [mesEmuList.size];

maxEmuSize = max(tmpModel.nCarbonMets(idNonPoolMets));
for s = 1 : maxEmuSize
    if any(mesEmuSizes==s)
        emuNetwork(s).emuList = mesEmuList(mesEmuSizes==s);
    else
        emuNetwork(s).emuList = [];        
    end
    emuNetwork(s).emuRxns = [];
end

%% �T�C�Y�����������Ȃ���EMU�l�b�g���[�N�����
for s = maxEmuSize:-1:1
    nTmpEmu = length(emuNetwork(s).emuList);
    if nTmpEmu > 0
        % �X�^�[�g��EMU���w�肷��
        for e = 1 : nTmpEmu
            tmpEmu = {emuNetwork(s).emuList(e).name};
            tmpMet = {emuNetwork(s).emuList(e).met};
            if tmpModel.massType(strcmp(tmpModel.mets, tmpMet)) <= 0
                emuNetwork = createEmuNetwork(tmpModel, tmpEmu, tmpMet, s, emuNetwork);
            end
        end
    end
end

tmpModel = rmfield(tmpModel, 'cellRxns');
tmpModel = rmfield(tmpModel, 'cellCTrans');
tmpModel = rmfield(tmpModel, 'arrowIds');


%% EMU�̐����J�E���g
nEmu = zeros(maxEmuSize,1);
nMassPattern = zeros(maxEmuSize,1);
for s = 1 : maxEmuSize
    if isempty(emuNetwork(s).emuList)
        continue
    end
    nEmu(s) = nnz([emuNetwork(s).emuList.isKnownEmu]==0);
    nMassPattern(s) = nEmu(s).*(s+1);
end
model.nEmu = nEmu;
model.nMassPattern = nMassPattern;
disp(['# EMU: ' num2str(sum(nEmu)) ', # mass pattern (or ODE): ' num2str(sum(nMassPattern))])
for s = 1 : maxEmuSize
    if nMassPattern(s) ==0
        continue
    end
    disp([' EMU size ' num2str(s) ': ' num2str(nMassPattern(s)) ' mass patterns (or ODEs)'])
end

model.emuNetwork= emuNetwork;

end

%% ���͑Ώۂ�EMU�ɘA������EMU���Ȃ�
function emuNetwork =  createEmuNetwork(model, emu, tmpMet, sizeEmu, emuNetwork)
% if sizeEmu == 3 && strcmp(emu, {'T3P123'})
%     keyboard;
% end

try
    metId = find(strcmp(model.mets, tmpMet));
    synsRxnIds = find(model.Sprod(metId,:)>0);
%     synsRxnIds = find(model.S(metId,:)>0);
    if isempty(synsRxnIds)
        return
    end
    prodEmu = emu;
    for sr = 1 : length(synsRxnIds)
        prodMet = tmpMet{:};
        prodMetCoef = model.Sprod(metId, synsRxnIds(sr));
%         prodMetCoef = model.S(metId, synsRxnIds(sr));
        % �O��̂�EMU��������
        [emu,tmpEmuList, idCorresMetInRxn, isMultiProd] = ...
            findPrecEmu(model, synsRxnIds(sr), prodMet, prodEmu, prodMetCoef);
        % EMU�̏����X�V
        [emuNetwork, isContFindEmu] = updateEmu(model, emuNetwork, tmpEmuList, ...
            sizeEmu, synsRxnIds(sr), prodEmu,  idCorresMetInRxn, isMultiProd);
        for j = 1 : length(emu)
            if isContFindEmu(j)
                % ����EMU��T��
                prodMet = emu{j}(1:end-sizeEmu);
                emuNetwork =  createEmuNetwork(model, emu(j), {prodMet}, sizeEmu, emuNetwork);
            end
        end
    end
catch err
%     keyboard;
    rethrow(err)
    return
end

end

%% �����̑O��̂ł���EMU��T��
function [emu, tmpEmuList, idCorresMetInRxn, isMultiProd] = ...
    findPrecEmu(model, synsRxnId, prodMet, prodEmu, prodMetCoef)

tmpCellRxn = model.cellRxns(synsRxnId,:);
tmpCellCTrans = model.cellCTrans(synsRxnId,:);
idMetInRxn = find(strcmp({prodMet}, tmpCellRxn));
idMetInRxn = idMetInRxn(idMetInRxn>model.arrowIds(synsRxnId)); % 5/26�X�V�Bproduct���̈ʒu�݂̂��C�ɂ���B
nProdMetInRxn = nnz(strcmp({prodMet}, tmpCellRxn(model.arrowIds(synsRxnId)+1:end)));

emu = cell(model.arrowIds(synsRxnId)/2,1);
idCorresMetInRxn = zeros(model.arrowIds(synsRxnId)/2,1);
emuMets = cell(model.arrowIds(synsRxnId)/2,1);
isCInEmu = cell(model.arrowIds(synsRxnId)/2,1);
subsEmuCoef =  zeros(model.arrowIds(synsRxnId)/2,1);
prodEmuCoef =  zeros(model.arrowIds(synsRxnId)/2,1);

jj = 0;
if length(idMetInRxn) >= 2
    isMultiProd = true;
else 
    isMultiProd = false;
end
for m = 1 : length(idMetInRxn)
    tmpCAtoms = [];
    % prodEmu�̒Y�f���q�ԍ� (abc...)�̎擾
    for j = 1 : length(prodEmu{:})-length(prodMet)
        idCarbon = str2double(prodEmu{1}(j+length(prodMet)));
        tmpCAtoms = [tmpCAtoms, tmpCellCTrans{idMetInRxn(m)}(idCarbon)];
    end
    % �O���Emu�̎擾    
    for j = 1 :model.arrowIds(synsRxnId)
%     for j = 1 :model.arrowIds(synsRxnId)/2
        tempMet = tmpCellRxn{j};
%         tempMet = tmpCellRxn{j*2-1};
        tmpIsCInEmu = ismember(tmpCellCTrans{j},  tmpCAtoms);
        if any(tmpIsCInEmu)
            jj = jj +1;
            isCInEmu{jj} = tmpIsCInEmu;            
            emu(jj) = defineEmu(tempMet, isCInEmu{jj});
            emuMets(jj) = {tempMet};
            idCorresMetInRxn(jj) = idMetInRxn(m);
            subsEmuCoef(jj) = 1;
            prodEmuCoef(jj) = prodMetCoef/nProdMetInRxn;
            if j >= 2 && strncmp(tmpCellRxn(j-1), {'('}, 1)
                subsEmuCoef(jj) = str2double(tmpCellRxn{j-1}(2:end-1));
            end
        end
    end
end
emuMets = emuMets(~cellfun('isempty',emu));
idCorresMetInRxn = nonzeros(idCorresMetInRxn);
isCInEmu = isCInEmu(~cellfun('isempty',emu));
emu = emu(~cellfun('isempty',emu));

for j = 1 : length(emu)
    tmpEmuList(j).name = emu{j};
    tmpEmuList(j).met = emuMets{j};
    tmpEmuList(j).isCInEmu = isCInEmu{j};
    tmpEmuList(j).size = nnz(isCInEmu{j});
    tmpEmuList(j).subsEmuCoef = subsEmuCoef(j);
    tmpEmuList(j).prodEmuCoef = prodEmuCoef(j);
end

% if prodEmuCoef ~=1
%     keyboard;
% end
end


%% EMU�̏����X�V����
function [emuNetwork, isContFindEmu] = updateEmu(...
    model, emuNetwork, tmpEmuList, sizeEmu, synsRxnId, prodEmu, idCorresMetInRxn, isMultiProd)

emu = {tmpEmuList.name};
emuMets = {tmpEmuList.met};

subsEmuCoef = [tmpEmuList.subsEmuCoef];
prodEmuCoef = [tmpEmuList.prodEmuCoef];
tmpEmuList = rmfield(tmpEmuList, {'subsEmuCoef', 'prodEmuCoef'});

isContFindEmu = false(length(emu),1);
isKnownEmu = true(length(emu),1);
for j = 1 : length(emu)
    tmpEmuList(j).isSmallerSize = false;
    tmpSizeEmu = tmpEmuList(j).size;
    
    % ��������EMU��ǉ�
    if ~ismember(emu(j), {emuNetwork(sizeEmu).emuList.name})
        if tmpSizeEmu == sizeEmu
            if model.massType(strcmp(model.mets, emuMets(j)))<=0
                isContFindEmu(j)=true;     % �or�v�[���̑�ӕ���EMU�Ȃ玟��T���Ȃ�
                isKnownEmu(j) =false;
            end
        else
            tmpEmuList(j).isSmallerSize = true;
        end
        tmpEmuList(j).isKnownEmu = isKnownEmu(j);
        emuNetwork(sizeEmu).emuList = [emuNetwork(sizeEmu).emuList, tmpEmuList(j)];
    end
    
    % �T�C�Y�̏�����EMU��ǉ�
    tmpEmuList(j).isKnownEmu  = true;
    tmpEmuList(j).isSmallerSize = false;
    if model.massType(strcmp(model.mets, emuMets(j)))<=0
        tmpEmuList(j).isKnownEmu  = false;
    end
    if isempty(emuNetwork(tmpSizeEmu).emuList)
        emuNetwork(tmpSizeEmu).emuList = [emuNetwork(tmpSizeEmu).emuList, tmpEmuList(j)];
    elseif ~ismember(emu(j), {emuNetwork(tmpSizeEmu).emuList.name})
        emuNetwork(tmpSizeEmu).emuList = [emuNetwork(tmpSizeEmu).emuList, tmpEmuList(j)];
    end
end

% ��������EMU���Ȃ��������L�^
if isMultiProd
% if length(emuMets)==2 &&  sum([tmpEmuList.size])~= sizeEmu %F6P �� T3P + T3P�̂悤�ȏꍇ�̗�O����
    [unqIdCorresMetInRxn, idUnq]= unique(idCorresMetInRxn);
    subsEmuCoef = subsEmuCoef(idUnq);
    prodEmuCoef = prodEmuCoef(idUnq);
    for i = 1 : length(unqIdCorresMetInRxn)
        tmpEmuRxns(i).subsEmu = emu(unqIdCorresMetInRxn(i)==idCorresMetInRxn);
        tmpEmuRxns(i).prodEmu = prodEmu;
        tmpEmuRxns(i).subsEmuCoef = subsEmuCoef(i);
        tmpEmuRxns(i).prodEmuCoef = prodEmuCoef(i);
        tmpEmuRxns(i).fluxId = synsRxnId;
    end
else
    tmpEmuRxns.subsEmu = emu;
    tmpEmuRxns.prodEmu = prodEmu;
    tmpEmuRxns.subsEmuCoef = subsEmuCoef;
    tmpEmuRxns.prodEmuCoef = prodEmuCoef(1);
    tmpEmuRxns.fluxId = synsRxnId;
end

if isempty(emuNetwork(sizeEmu).emuRxns)
    emuNetwork(sizeEmu).emuRxns = tmpEmuRxns;
else
    emuNetwork(sizeEmu).emuRxns = [emuNetwork(sizeEmu).emuRxns, tmpEmuRxns];
end
    
end


%% EMU���`����
function emu = defineEmu(met, isCInEmu)
idEmu = 0;
for k = length(isCInEmu) : -1 : 1
    if isCInEmu(k)
        idEmu = idEmu*10+k;
    end
end
idEmuStr = num2str(idEmu);
idEmuStr = idEmuStr(end:-1:1);
emu = {[met idEmuStr]};
end

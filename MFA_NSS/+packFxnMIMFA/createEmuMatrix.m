function model = createEmuMatrix161129(model, optionsMFA)
%createEmuMatrix EMU�̎��ԕω������߂�s����쐬
% dX/dt = {(A - S'�Ev)X + B�EY}/C 
% X: ���mEMU
% Y: ���mEMU
% S': ���w�ʘ_�s��i�������AEMU�̑�ӕ��Ɋւ�����̂̂݁j
% A, B: �t���b�N�X�Ɉˑ����Č��܂�s��
%   A = A'1*v1 + A'2*v2 + A'3*v3 + ...
%   A = B'1*v1 + B'2*v2 + B'3*v3 + ...

[~, nRxns] = size(model.S);
emuNetwork = model.emuNetwork;

for s = 1 : length(emuNetwork)
    emuList = emuNetwork(s).emuList;
    emuRxns = emuNetwork(s).emuRxns;
    if isempty(emuList)
        continue
    end
%     if isempty(emuRxns)
%         continue
%     end
    
    %% ������Emu���獇�������ꍇ�ɂ��ď���
    synsEmuList = [];
    nSynsEmu = 0;
    if isempty(emuRxns)
        prodEmu = [];
        subsEmu = [];
        emu = {emuList([emuList.size]==s).name}; %Pyr+CO2-->OAA�ǉ����̃o�O����
    else
        for j = 1 : length(emuRxns)
            if length(emuRxns(j).subsEmu) > 1
                nSynsEmu = nSynsEmu+1;
                synsEmuName = ['synsEmu' num2str(nSynsEmu) '_' emuRxns(j).subsEmu{:}];
                synsEmuList(nSynsEmu).name = synsEmuName;
                synsEmuList(nSynsEmu).correctEmu = emuRxns(j).subsEmu;
                emuRxns(j).subsEmu = {synsEmuName};
                
                [~, idEmu] = ismember(synsEmuList(nSynsEmu).correctEmu, {emuList.name});
                synsEmuList(nSynsEmu).correctEmuSize = [emuList(idEmu).size];
            end
        end
        
        prodEmu = [emuRxns.prodEmu];
        subsEmu = [emuRxns.subsEmu];
        emu = union({emuList([emuList.size]==s).name},unique([subsEmu, prodEmu])); %Pyr+CO2-->OAA�ǉ����̃o�O����
    end
%     emu1 = union({emuList([emuList.size]==s).name},unique([subsEmu, prodEmu]));
%     emu2 = unique([subsEmu, prodEmu]);
    
    %% ���m��Emu�����m��Emu�̏��ɕ��בւ�
    isKnownEmu = [emuList.isKnownEmu];
    unknownEmu = intersect(emu,{emuList(~isKnownEmu).name});
    knownEmu = setdiff(emu, unknownEmu);
    emu = [unknownEmu, knownEmu];
    nUnkEmu = length(unknownEmu);
    nKnownEmu = length(knownEmu);

    isNotSingleEmu = false(1,nKnownEmu);
    knownEmuSize = cell(nKnownEmu, 1);
    if ~ isempty(synsEmuList)
        synsEmu = {synsEmuList.name};
        [isNotSingleEmu, idIsMember] = ismember(knownEmu, synsEmu);
        knownEmu(isNotSingleEmu) = {synsEmuList(idIsMember(isNotSingleEmu)).correctEmu};
        knownEmuSize (isNotSingleEmu) =  {synsEmuList(idIsMember(isNotSingleEmu)).correctEmuSize};
    end
    isSingleEmu = ~isNotSingleEmu';
    
    %% A��B�̍쐬
    for j =  1 : nRxns
        if isempty(emuRxns)
            emuNetwork(s).arrayMakeA.array{j} =[];
            emuNetwork(s).arrayMakeB.array{j} = [];
        else
            tmpArray = zeros(nUnkEmu, length(emu));
            idEmuRxns = find([emuRxns.fluxId] == j);
            if ~ isempty(idEmuRxns)
                prodEmus =  [emuRxns(idEmuRxns).prodEmu];
                subsEmus = [emuRxns(idEmuRxns).subsEmu];
                %             idr = find(ismember(emu, prodEmus));
                %             idc =  find(ismember(emu, subsEmus));
                for k = 1 : length(prodEmus)
                    idr = find(strcmp(prodEmus(k), emu));
                    idc = find(strcmp(subsEmus(k), emu));
                    tmpArray(idr,idc) = tmpArray(idr,idc) +emuRxns(idEmuRxns(k)).prodEmuCoef;
                    %                 tmpArray(idr,idc) = tmpArray(idr,idc) +1;
                end
            end
            idConsMet = find(model.Ssubs(:,j)<0);
            coefs = model.Ssubs(idConsMet,j);
            for l = 1 : nnz(idConsMet)
                isConsEmu = ismember({emuList.met}, model.mets(idConsMet(l)));
                isConsUnkEmu = ismember(unknownEmu, {emuList(isConsEmu).name});
                if any(isConsUnkEmu)
                    id = find(isConsUnkEmu);
                    for k = 1 : length(id)
                        tmpArray(id(k),id(k)) = tmpArray(id(k),id(k))+coefs(l);  %�����炭����ł悢�̂ł͂Ȃ��ł��傤���H
                    end
                end
            end
            
            emuNetwork(s).arrayMakeA.array{j} = sparse(tmpArray(:,1:nUnkEmu));
            emuNetwork(s).arrayMakeB.array{j} = sparse(tmpArray(:,nUnkEmu+1:end));
        end
    end
    
    % metabolite ID of EMU
    idMetA = zeros(nUnkEmu, 1);
    for e = 1 : nUnkEmu
        idEmuList = strcmp(unknownEmu(e), {emuList.name});
        idMetA(e) = find(strcmp({emuList(idEmuList).met}, model.mets));
    end
    idMetB = zeros(nKnownEmu, 1);
    for e = 1 : nKnownEmu
        idEmuList = strcmp(knownEmu(e), {emuList.name});
        if any(idEmuList)
            idMetB(e) = find(strcmp({emuList(idEmuList).met}, model.mets));
        end
    end
    
    emuNetwork(s).arrayMakeA.emu = unknownEmu';
    emuNetwork(s).arrayMakeA.idMet = idMetA;
    emuNetwork(s).arrayMakeB.emu = knownEmu';
    emuNetwork(s).arrayMakeB.idMet = idMetB;
    emuNetwork(s).arrayMakeB.isSingleEmu = isSingleEmu;
    emuNetwork(s).arrayMakeB.sizeEmu = knownEmuSize;
    
    %% S�̂����AEmu���R�������ӕ��ɂ��Ă̂ݔ����o���B
    [~,idEmuList] = ismember(unknownEmu, {emuList.name});
    [~, idMet] = ismember({emuList(idEmuList).met}, model.mets);
    emuNetwork(s).emuS = model.S(idMet,:);
    
    %% ��������A, B���쐬 (�m�F�p)
    cellArrayA = cell(nUnkEmu);
    cellArrayB = cell(nUnkEmu,nKnownEmu);
    for j =  1 : nRxns
        id = find(emuNetwork(s).arrayMakeA.array{j}~=0);
        for k = 1 : length(id)
            preCoef = emuNetwork(s).arrayMakeA.array{j}(id(k));
            switch preCoef
                case 1
                    coef = ' + ';
                case -1
                    coef = ' - ';
                otherwise
                    if preCoef > 0
                        coef = [' + ' num2str(preCoef) ' '];
                    else
                        coef = [' ' num2str(preCoef) ' '];
                    end
            end
            cellArrayA{id(k)} = [cellArrayA{id(k)} coef  'v' num2str(j)];
        end
        
        id = find(emuNetwork(s).arrayMakeB.array{j}~=0);
        for k = 1 : length(id)
            preCoef = emuNetwork(s).arrayMakeB.array{j}(id(k));
            switch preCoef
                case 1
                    coef = ' + ';
                case -1
                    coef = ' - ';
                otherwise
                    if preCoef > 0
                        coef = [' + ' num2str(preCoef) ' '];
                    else
                        coef = [' - ' num2str(preCoef) ' '];
                    end
            end
            cellArrayB{id(k)} = [cellArrayB{id(k)} coef  'v' num2str(j)];
        end
    end
%     emuNetwork(i).arrayMakeA.cellArrayA=cellArrayA;
%     emuNetwork(i).arrayMakeB.cellArrayB=cellArrayB;

model.emuNetwork= emuNetwork;
end


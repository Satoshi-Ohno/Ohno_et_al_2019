%% エクセルデータより情報整理
function model = loadXlsNetwork(model,optionsMFA, fileName)
%loadXlsNetwork  Load metabolic network defined in Excel file
% 
% model = loadXlsNetwork(fileName)
% 
% 

%% Load reaction information
[~,~,rawXlsData] = xlsread(fileName,'Reactions');
% [~,~,rawXlsData] = xlsread(fileName,'MassBalance');

colNames = rawXlsData(1,:);
rawXlsData = rawXlsData(2:end,:);
for i = 1 : size(rawXlsData,2);
    switch colNames{i}
        case {'Reaction IDs'}
            tmpRxnInfo.fluxIDs = rawXlsData(:,i);
        case {'Reaction abbreviations'}
            tmpRxnInfo.rxnNames = rawXlsData(:,i);
        case {'Reaction names'}
            tmpRxnInfo.rxnFullNames = rawXlsData(:,i);
        case {'Reversible reaction'}
            tmpRxnInfo.isRev = cell2mat(rawXlsData(:,i));
        case {'Reaction formula'}
            tmpRxnInfo.rxns = rawXlsData(:,i);
        case {'Carbon atom transitions'}
            tmpRxnInfo.carbonTransitions = rawXlsData(:,i);
        case {'Note'}
            tmpRxnInfo.note = rawXlsData(:,i);
    end
end
nTmpRxns = length(tmpRxnInfo.rxnNames);

% fieldNames = rawXlsData(1,:);
% rawXlsData = rawXlsData(2:end,:);
% for i = 1 : size(rawXlsData,2);
%     switch fieldNames{i}
%         case {'revSets', 'minInitFlux', 'maxInitFlux'}
%             model.rxnInfo.(fieldNames{i}) = cell2mat(rawXlsData(:,i));
%         otherwise
%         model.rxnInfo.(fieldNames{i}) = rawXlsData(:,i);
%     end
% end

nRxns = nTmpRxns + nnz(tmpRxnInfo.isRev);
rxnInfo.fluxIds = zeros(nRxns,1);
rxnInfo.rxnNames = cell(nRxns,1);
rxnInfo.rxnFullNames = cell(nRxns,1);
rxnInfo.revSets = zeros(nRxns,1);
rxnInfo.rxns = cell(nRxns,1);
rxnInfo.carbonTransitions = cell(nRxns,1);
rxnInfo.note = cell(nRxns,1);

ii = 0;
idRevSets = 0;
for i = 1 : nTmpRxns
    if tmpRxnInfo.isRev(i)
        idRevSets = idRevSets + 1;        
        % forward
        ii = ii + 1;
        rxnInfo.fluxIds(ii) = ii;
        rxnInfo.rxnNames{ii} =  [tmpRxnInfo.rxnNames{i} '_f'];
        rxnInfo.rxnFullNames{ii} =  [tmpRxnInfo.rxnFullNames{i} ' (forward)'];
        rxnInfo.revSets(ii) = idRevSets;
        rxnInfo.rxns(ii) =  tmpRxnInfo.rxns(i);
        rxnInfo.carbonTransitions(ii) = tmpRxnInfo.carbonTransitions(i);
        rxnInfo.note(ii) = tmpRxnInfo.note(i);
        
        % backward
        ii = ii + 1;
        rxnInfo.fluxIds(ii) = ii;
        rxnInfo.rxnNames{ii} =  [tmpRxnInfo.rxnNames{i} '_r'];
        rxnInfo.rxnFullNames{ii} =  [tmpRxnInfo.rxnFullNames{i} ' (backward)'];
        rxnInfo.revSets(ii) = idRevSets;
        tmp = strsplit(tmpRxnInfo.rxns{i});
        tmp = tmp(end:-1:1);
        rxnInfo.rxns{ii} = strjoin(tmp, ' ');
%         rxnInfo.rxns{ii} = [tmp{:}];
        tmp = strsplit(tmpRxnInfo.carbonTransitions{i});
        tmp = tmp(end:-1:1);
        rxnInfo.carbonTransitions{ii} =strjoin(tmp, ' ');
%         rxnInfo.carbonTransitions{ii} = [tmp{:}];
        rxnInfo.note(ii) = tmpRxnInfo.note(i);
    else
        ii = ii + 1;
        rxnInfo.fluxIds(ii) = ii;
        rxnInfo.rxnNames(ii) =  tmpRxnInfo.rxnNames(i);
        rxnInfo.rxnFullNames(ii) =  tmpRxnInfo.rxnFullNames(i);
        rxnInfo.revSets(ii) = 0;
        rxnInfo.rxns(ii) =  tmpRxnInfo.rxns(i);
        rxnInfo.carbonTransitions(ii) = tmpRxnInfo.carbonTransitions(i);
        rxnInfo.note(ii) = tmpRxnInfo.note(i);        
    end
    
end
model.rxnInfo = rxnInfo;

%% Load metabolite information
[~,~,xlsMetData] = xlsread(fileName,'Metabolites');
% [N,C] = xlsread(fileName,'Metabolites');
% xlsMetData = numIcell(N, C, [2,3]);
colNames = xlsMetData(1,:);
% fieldNames = xlsMetData(1,:);
for i = 1 : size(xlsMetData,2)
    switch colNames{i}
        case {'Metabolite IDs'}
            model.metInfo.metIds = xlsMetData(2:end, i);
        case {'Metabolite abbreviations'}
            model.metInfo.mets =  xlsMetData(2:end, i);
        case {'Metabolite names'}
            model.metInfo.metFullNames =  xlsMetData(2:end, i);
        case {'KEGG IDs'}
            model.metInfo.metKeggIds = xlsMetData(2:end, i);
        case {'# of skeletal carbon atoms'}
            model.metInfo.nCarbonMets = cell2mat(xlsMetData(2:end, i));
        case {'Metabolite concentrations used for calculating wRSS'}
            model.metInfo.isEvalConc = logical(cell2mat(xlsMetData(2:end, i)));
        case {'Mass isotopomer fractions used for calculating wRSS'}
            model.metInfo.massType = ~logical(cell2mat(xlsMetData(2:end, i)));
        case {'Outside the metabolic network'}
            model.metInfo.isPoolMets = logical(cell2mat(xlsMetData(2:end, i)));
        case {'Note'}
            model.metInfo.note = xlsMetData(2:end, i);
    end
end

% model.isPoolMets = logical(cell2mat(model.isPoolMets));
% model.compt = cell2mat(model.compt);
% model.isSameEval = cell2mat(model.isSameEval);
% model.isEvalConc = logical(cell2mat(model.isEvalConc));
% model.massType = cell2mat(model.massType);
% model.nCarbonMets = cell2mat(model.nCarbonMets);

model.rxnNames = model.rxnInfo.rxnNames;
model.rxns = model.rxnInfo.rxns;
model.mets = model.metInfo.mets;
nMets = length(model.mets);


model.metInfo.compt = zeros(nMets,1)+1;
model.metInfo.idSameEval = zeros(nMets,1);

model.metInfo.massType = double(model.metInfo.massType);
model.metInfo.massType(model.metInfo.massType==1) = 2;
idGlc_ex = find(strcmp({'Glc_ex'}, model.metInfo.mets));
model.metInfo.massType(idGlc_ex) = 1;


%% Sperate merged kegg IDs
metKeggIds = cell(nMets, 3);
for i = 1 : nMets
    if isnan(model.metInfo.metKeggIds{i})
        continue
    end
    tmpKeggId = model.metInfo.metKeggIds{i};
    splitKeggId = strsplit(tmpKeggId, '_');
    metKeggIds(i,1:length(splitKeggId)) = splitKeggId;
end
model.metInfo.metKeggIds = metKeggIds;

%% Load metabolite formula from KEGG
base = 'http://rest.kegg.jp/';
operation = 'get/';
model.metInfo.formulas = cell(size(model.metInfo.metKeggIds));
for i = 1 : numel(model.metInfo.metKeggIds)
    if isempty(model.metInfo.metKeggIds{i})
        continue
    end
    cpdInfo = urlread(strcat(base,operation,model.metInfo.metKeggIds{i}));
    locStartFormula = regexpi(cpdInfo, 'FORMULA');
    locEndFormula = regexpi(cpdInfo, 'EXACT_MASS')-2;
    model.metInfo.formulas{i} = cpdInfo(locStartFormula+12:locEndFormula);
end

Ematrix = cell(size(model.metInfo.metKeggIds));
model.metInfo.atoms = {'C','N', 'O', 'H', 'P', 'Others'};
for i = 1 : numel(model.metInfo.metKeggIds)
    if isempty(model.metInfo.metKeggIds{i})
        continue
    end
    tmpEmatrix = zeros(1,length(model.metInfo.atoms)-1);
    for j = 1 : length(model.metInfo.atoms)-1
        locStartAtom = find(model.metInfo.atoms{j}==model.metInfo.formulas{i});
        if isempty(locStartAtom)
            continue
        end
        isAtom = isletter(model.metInfo.formulas{i});
        locEndAtom =find(isAtom(locStartAtom+1:end),1);
        if isempty(locEndAtom)
            if locStartAtom == length(model.metInfo.formulas{i})
                tmpEmatrix(j) = 1;
            else
                tmpEmatrix(j) = str2double(model.metInfo.formulas{i}(locStartAtom+1:end));
            end
        elseif locEndAtom == 1
            tmpEmatrix(j) = 1;
        else
            tmpEmatrix(j) = str2double(model.metInfo.formulas{i}(locStartAtom+1:locStartAtom+locEndAtom-1));
        end
    end
    Ematrix{i} = tmpEmatrix;
end
model.metInfo.Ematrix=Ematrix;

%% Make a matrix to represent mass isotopomer fractions used for calculating wRSS
useMassMat = false(nMets, max(model.metInfo.nCarbonMets)+1);
for i = 1 : nMets
    if model.metInfo.massType(i)
        continue
    end
    useMassMat(i,1:model.metInfo.nCarbonMets(i)) = true;
end
%     keyboard;
%     
%     if isempty(model.metInfo.useMass{i})
%         continue
%     end
%     if isnan(model.metInfo.useMass{i})
%         continue
%     end
%     useMass = strsplit(model.metInfo.useMass{i}, '_');
%     for j = 1:length(useMass)
%         tmpMass = str2double(useMass{j}(2:end));
%         useMassMat(i,tmpMass+1) = true;
%     end
% end
model.metInfo.useMassMat = useMassMat;

%% out of calibration rangeのデータをロード
model.outCalib.idMets = [];
model.outCalib.expIds= [];
model.outCalib.mass = [];
model.outCalib.time = [];

%%  基質の質量同位体のうちの同位体割合の定義
% [~,~,rawCell] = xlsread(fileName,'CarbonSource');
% rawCell = rawCell(2:end,:);
% for n = 1 : numel(rawCell)
%     if isnan(rawCell{n})
%         rawCell{n} = [];
%     end
% end
% 
% [~,idUnqMet,idRawMet] = unique(rawCell(:,2));
% nSubs = length(idUnqMet);
% subsIsoRatio.idMets=zeros(nSubs,1);
% subsIsoRatio.mets=cell(nSubs,1);
% subsIsoRatio.mass=cell(nSubs,1);
% subsIsoRatio.CAtom13C=cell(nSubs,1);
% subsIsoRatio.isoRatio=cell(nSubs,1);
% subsIsoRatio.isoMat=cell(nSubs,1);
% subsIsoRatio.matMDV2IDV=cell(nSubs,1);
% 
% % 各代謝物についてまとめる。
% for m = 1 : length(idUnqMet)    
%     tmpRawCell = rawCell(idRawMet==m,:);
%     idMet = find(strcmp(tmpRawCell(1,2), model.mets));
%     if isempty(idMet)
%         error('metabolite name of the sheet "SubstrateIsotopomerRatio" is not correct.')
%     end
%     subsIsoRatio.idMets(m) = idMet;
%     subsIsoRatio.mets{m} = model.mets(idMet);
%     
%     % 各質量代謝物についてまとめる
%     [~, idUnqMass, idRawMass]= unique(tmpRawCell(:,3));    
%     nTmpMDV = length(idUnqMass);
%     subsIsoRatio.mass{m} = zeros(nTmpMDV,1);
%     subsIsoRatio.CAtom13C{m} = cell(nTmpMDV,1);
%     subsIsoRatio.isoRatio{m} = cell(nTmpMDV,1);
%     for mi = 1 : nTmpMDV
%         loc = idRawMass == mi;
%         subsIsoRatio.mass{m}(mi) = str2double(tmpRawCell{idUnqMass(mi),3}(2:end));
%         subsIsoRatio.CAtom13C{m}{mi} = tmpRawCell(loc,4);
%         subsIsoRatio.isoRatio{m}{mi} = tmpRawCell(loc,5);
%     end
%     mass = subsIsoRatio.mass{m};
%     
%     % 今回考慮するすべての同位体の行列を作る。
%     % 行：同位体、列：炭素原子
%     % 1のところが13C
%     CAtom13CList = vertcat(subsIsoRatio.CAtom13C{m}{:});
%     nIDV = length(CAtom13CList);
%     CAtomList = unique([CAtom13CList{:}]);
%     nCarbonAtom = length(CAtomList);    
%     isoMat = false(nIDV, nCarbonAtom);
%     for i = 1 : nIDV;
%         if isempty(CAtom13CList{i})
%             continue
%         end
%         isoMat(i,:) = ismember(CAtomList, CAtom13CList{i});
%     end
%     subsIsoRatio.isoMat{m} = isoMat;
%     
%     % 質量同位体割合(MDV)から同位体割合 (IDV)を作る行列の作成
%     nMDV = nCarbonAtom+1;
%     subsIsoRatio.matMDV2IDV{m} = zeros(nIDV,nMDV);
%     for mi = 1 : nTmpMDV
%         locIDV =  sum(isoMat,2)==mass(mi);
%         tmpRatio = subsIsoRatio.isoRatio{m}{mi};
%         subsIsoRatio.matMDV2IDV{m}(locIDV,mass(mi)+1) = cell2mat(tmpRatio);
%     end
% end
% 
% model.metInfo.subsIsoRatio = subsIsoRatio;

end





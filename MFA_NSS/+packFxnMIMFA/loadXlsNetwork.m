%% Load network defined in Excel file
function model = loadXlsNetwork(model,optionsMFA, fileName)

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
        tmp = strsplit(tmpRxnInfo.carbonTransitions{i});
        tmp = tmp(end:-1:1);
        rxnInfo.carbonTransitions{ii} =strjoin(tmp, ' ');
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
colNames = xlsMetData(1,:);
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

model.rxnNames = model.rxnInfo.rxnNames;
model.rxns = model.rxnInfo.rxns;
model.mets = model.metInfo.mets;
nMets = length(model.mets);


model.metInfo.idSameEval = zeros(nMets,1);

% massType
% 0: simulated
% 1: always U-13C
% 2: always U-12C
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
model.metInfo.useMassMat = useMassMat;

%% out of calibration rangeのデータをロード
model.outCalib.idMets = [];
model.outCalib.expIds= [];
model.outCalib.mass = [];
model.outCalib.time = [];

end





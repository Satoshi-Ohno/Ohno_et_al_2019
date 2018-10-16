%% �����f�[�^�̑�ӕ��Ƃ̑Ή����쐬
function expData = connectExpModelMets170105(expData, preExpData, model)

%% �����f�[�^�Ƃ�KEGG ID�őΉ�

expData.mets = model.mets;
expData.nCarbonMets = model.metInfo.nCarbonMets;

expData.metKeggIds = model.metInfo.metKeggIds;
% expMets = [preExpData.cpdName]';
expMets = {preExpData.cpdName}';
if ischar(preExpData(1).keggId)
    preExpMetKeggIds = {preExpData.keggId}';
else
    preExpMetKeggIds = [preExpData.keggId]';
end
expData.expMets = cell(size(model.metInfo.metKeggIds));
for i = 1 : size(model.metInfo.metKeggIds,2)
    idNotEmp = find(~cellfun('isempty', model.metInfo.metKeggIds(:,i)));
    [isMesMets, idIsMember] = ismember(model.metInfo.metKeggIds(idNotEmp,i), preExpMetKeggIds);
    expData.expMets(idNotEmp(isMesMets),i) = expMets(idIsMember(isMesMets));
end


end
function vecExpConcStruct = makeVecExpConcs(model, expData, optionsMFA)

field2var(optionsMFA.varSet);

%% �����f�[�^�ƕ��U�̃x�N�g��
nExpConcs  = nEvalConcMets*nExpTime;
vecExpConcs = zeros(nExpConcs,1);
vecSDs = zeros(nExpConcs,1);
vecIdMet = zeros(nExpConcs,1);
vecIdTime = zeros(nExpConcs,1);
vecIsComp = zeros(nExpConcs,1);
mm = 0;
for m = 1 : nEvalConcMets
    idMet = idEvalConcMets(m);
        
    %��r����Conc�̑I��
    isComp = expData.concs(idMet,:)>0;  % 0����thresholdMDV�ɏC��
%     nComp = nnz(isComp);

    % �����̑�ӕ���merge���ĕ]������ꍇ(AcCoA��cyto��mito�Ȃ�)
    if idSameEval(idMet)>=1 && idMet >=2 && any(idSameEval(1:idMet-1)==idSameEval(idMet))
        isComp(:) = false;
    end
    
    % �����̑�ӕ��Z�x��SD�f�[�^�̃x�N�g�����쐬
    tmpExpConcs = reshape(expData.concs(idMet, :), nExpTime, 1);
    tmpSDs = reshape(expData.concsSD(idMet, :), nExpTime, 1);
%     tmpExpConcs = reshape(expData.concs(idMet, isComp), nComp, 1);
%     tmpSDs = reshape(expData.concsSD(idMet, isComp), nComp, 1);
    
    %�����_���̏C��
    tmpExpConcs(isnan(tmpExpConcs))=0;
    tmpSDs(tmpExpConcs==0)=0;
    
%     % compartment�ɂ��␳������͂��ƂŎ�������B
%     [tmpExpConcs, tmpSDs] = corrComptMedia2Cell(model, expData, optionsMFA, ...
%         tmpExpConcs, tmpSDs, coefCorrCompt, idMet);

    % �܂Ƃ�
    vecExpConcs(mm+(1:nExpTime)) = tmpExpConcs;
    vecSDs(mm+(1:nExpTime)) = tmpSDs;
    vecIdMet(mm+(1:nExpTime)) = idMet;
    vecIdTime(mm+(1:nExpTime)) = (1:nExpTime)';
    vecIsComp(mm+(1:nExpTime)) = isComp;
    
    mm = mm + nExpTime;
end
nExpConcs = mm;

vecExpConcStruct.vecExpConcs = vecExpConcs;
vecExpConcStruct.vecSDs = vecSDs;
vecExpConcStruct.vecIdMet = vecIdMet;
vecExpConcStruct.vecIdTime = vecIdTime;
vecExpConcStruct.vecIsComp = vecIsComp;
vecExpConcStruct.nExpConcs = nExpConcs;

end
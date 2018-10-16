function vecExpConcStruct = makeVecExpConcs(model, expData, optionsMFA)

field2var(optionsMFA.varSet);

%% 実験データと分散のベクトル
nExpConcs  = nEvalConcMets*nExpTime;
vecExpConcs = zeros(nExpConcs,1);
vecSDs = zeros(nExpConcs,1);
vecIdMet = zeros(nExpConcs,1);
vecIdTime = zeros(nExpConcs,1);
vecIsComp = zeros(nExpConcs,1);
mm = 0;
for m = 1 : nEvalConcMets
    idMet = idEvalConcMets(m);
        
    %比較するConcの選択
    isComp = expData.concs(idMet,:)>0;  % 0からthresholdMDVに修正
%     nComp = nnz(isComp);

    % 複数の代謝物をmergeして評価する場合(AcCoAのcytoとmitoなど)
    if idSameEval(idMet)>=1 && idMet >=2 && any(idSameEval(1:idMet-1)==idSameEval(idMet))
        isComp(:) = false;
    end
    
    % 実験の代謝物濃度とSDデータのベクトルを作成
    tmpExpConcs = reshape(expData.concs(idMet, :), nExpTime, 1);
    tmpSDs = reshape(expData.concsSD(idMet, :), nExpTime, 1);
%     tmpExpConcs = reshape(expData.concs(idMet, isComp), nComp, 1);
%     tmpSDs = reshape(expData.concsSD(idMet, isComp), nComp, 1);
    
    %欠損点等の修正
    tmpExpConcs(isnan(tmpExpConcs))=0;
    tmpSDs(tmpExpConcs==0)=0;
    
%     % compartmentによる補正→これはあとで実装する。
%     [tmpExpConcs, tmpSDs] = corrComptMedia2Cell(model, expData, optionsMFA, ...
%         tmpExpConcs, tmpSDs, coefCorrCompt, idMet);

    % まとめ
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
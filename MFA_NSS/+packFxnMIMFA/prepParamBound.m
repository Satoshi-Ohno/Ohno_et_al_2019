function [lb, ub] = prepParamBound170913(model, expData, optionsMFA, optimType, expandRatio)

if nargin <= 4
    expandRatio = 1;
end

field2var(optionsMFA.varSet)

if ~optionsMFA.isUseQPInMH  % 180830�C��
    if strcmp(optimType, 'metaheuristic')
        optimType = 'local';
    end
end

switch optimType
    case {'metaheuristic'}
        nOptimMetsInitConcs =  nOuterOptimMetsInitConcs;
        idOptimMetsInitConcs =idOuterOptimMetsInitConcs;
        nOptimMetsConcRates = nOuterOptimMetsConcRates;
        idOptimMetsConcRates = idOuterOptimMetsConcRates;
    case {'init', 'local'}
        nOptimMetsInitConcs = nNonPoolMets;
        idOptimMetsInitConcs =idNonPoolMets;
        nOptimMetsConcRates = nNonPoolMets;
        idOptimMetsConcRates =idNonPoolMets;
end


%% initConcs, knotConcRates, knotFluxes, knots��bound
lb.initConcs= zeros(nOptimMetsInitConcs, 1)+optionsMFA.lb.initConcs;
lb.knotConcRates= zeros(nOptimMetsConcRates, nKnots+2)+optionsMFA.lb.knotConcRates;
lb.knotFluxes = zeros(nIndFluxes, nKnots+2)+optionsMFA.lb.knotFluxes;
lb.knots =  zeros(nKnots,1)+optionsMFA.lb.knots+optionsMFA.minKnotTimeDiff;

ub.initConcs= zeros(nOptimMetsInitConcs, 1)+optionsMFA.ub.initConcs;
ub.knotConcRates= zeros(nOptimMetsConcRates, nKnots+2)+optionsMFA.ub.knotConcRates;
ub.knotFluxes = zeros(nIndFluxes, nKnots+2)+optionsMFA.ub.knotFluxes;
ub.knots =  zeros(nKnots,1)+optionsMFA.ub.knots-optionsMFA.minKnotTimeDiff;

switch optimType
    case {'metaheuristic'}
        if optionsMFA.isIncludeAllInitConcsInParamMH
            loc = isEvalConcMets(1:nNonPoolMets);
            lb.initConcs(loc,:) = optionsMFA.lb.mesInitConcs;
            ub.initConcs(loc,:) = optionsMFA.ub.mesInitConcs;
%             loc = ismember(model.compt, optionsMFA.idComptInnerOptim) & isEvalConcMets;
%             loc = loc(1:nNonPoolMets);
%             lb.initConcs(loc,:) = optionsMFA.lb.mesInitConcs;
%             ub.initConcs(loc,:) = optionsMFA.ub.mesInitConcs;
        end
    case {'init', 'local'}
        % �����ӕ� (MH�ł�QP�t�B�b�g�Alocal�Ŕ���`�t�B�b�g)��bound
        loc = isEvalConcMets(1:nNonPoolMets);
        lb.initConcs(loc,:) = optionsMFA.lb.mesInitConcs;
        lb.knotConcRates(loc,:) = optionsMFA.lb.mesKnotConcRates;
        ub.initConcs(loc,:) = optionsMFA.ub.mesInitConcs;
        ub.knotConcRates(loc,:) = optionsMFA.ub.mesKnotConcRates;
end


% ������ӕ��Z�x�ɂ��Ă̐������
if optionsMFA.isInitConcConstr
    tmpConstr = optionsMFA.constr.initConcs;
    [~, idTmp] = ismember(tmpConstr.mets, model.mets);
    for i = 1 : length(idTmp)
        lb.initConcs(idOptimMetsInitConcs==idTmp(i)) = tmpConstr.lb(i);
        ub.initConcs(idOptimMetsInitConcs==idTmp(i)) = tmpConstr.ub(i);
    end
end
% ��ӕ��Z�x�̕ω����x�ɂ��Ă̐������
if optionsMFA.isConcRateConstr
    tmpConstr = optionsMFA.constr.concRates;
    [~, idTmp] = ismember(tmpConstr.mets, model.mets);
    for i = 1 : length(idTmp)
        lb.knotConcRates(idOptimMetsConcRates==idTmp(i),:) = tmpConstr.lb(i);
        ub.knotConcRates(idOptimMetsConcRates==idTmp(i),:) = tmpConstr.ub(i);
    end
end

%% ���ʂ���bound
% coefCorrCompt
if isempty(idComptCorrParam)
    lb.coefCorrCompt = [];
else
    lb.coefCorrCompt = zeros(nComptCorrParam,1)+optionsMFA.lb.coefCorrCompt;
end
if isempty(idComptCorrParam)
    ub.coefCorrCompt = [];
else
    ub.coefCorrCompt = zeros(nComptCorrParam,1)+optionsMFA.ub.coefCorrCompt;
end

%�����l���̏���
if optionsMFA.isInitSteadyState
    ub.knotConcRates(:,1) = +10^-6;
    lb.knotConcRates(:,1) = -10^-6;
end

% �t���b�N�X�ɂ��Ă̐������
if optionsMFA.isFluxConstr
    tmpConstr = optionsMFA.constr.fluxes;
    [~, idTmp] = ismember(tmpConstr.idRxns, model.idIndFluxes);
    for i = 1 : length(idTmp)
        if idTmp(i)==0
            continue
        end
        lb.knotFluxes(idTmp(i),tmpConstr.idFullKnot(i)) = tmpConstr.lb(i);
        ub.knotFluxes(idTmp(i),tmpConstr.idFullKnot(i)) = tmpConstr.ub(i);
    end
end


%% expandRatio�ɏ]���ăp�����[�^�͈͂��g��
lb.initConcs = lb.initConcs/expandRatio;
lb.knotConcRates= lb.knotConcRates*expandRatio;
lb.knotFluxes = lb.knotFluxes/expandRatio;
lb.knots =  lb.knots -optionsMFA.minKnotTimeDiff+optionsMFA.minKnotTimeDiff/expandRatio;

ub.initConcs= ub.initConcs*expandRatio;
ub.knotConcRates= ub.knotConcRates*expandRatio;
ub.knotFluxes = ub.knotFluxes*expandRatio;
ub.knots =  ub.knots +optionsMFA.minKnotTimeDiff-optionsMFA.minKnotTimeDiff/expandRatio;

if ~isempty(idComptCorrParam)
    if  lb.coefCorrCompt/expandRatio <  ub.coefCorrCompt*expandRatio
        lb.coefCorrCompt = lb.coefCorrCompt/expandRatio;
        ub.coefCorrCompt = ub.coefCorrCompt*expandRatio;
    else
        tmpCoefCorrCompt = (lb.coefCorrCompt + ub.coefCorrCompt)/2;
        lb.coefCorrCompt = tmpCoefCorrCompt + expandRatio*(lb.coefCorrCompt - tmpCoefCorrCompt);
        ub.coefCorrCompt = tmpCoefCorrCompt + expandRatio*(ub.coefCorrCompt - tmpCoefCorrCompt);
    end
end

%% �S�t���b�N�X�ɂ��Ă̐������
switch optimType
    case {'init','local'}
        lb.knotFluxesAll = zeros(nRxns,nKnots+2)+optionsMFA.lb.knotFluxes;
        ub.knotFluxesAll = zeros(nRxns,nKnots+2)+optionsMFA.ub.knotFluxes;
        % �t���b�N�X�ɂ��Ă̐������
        if optionsMFA.isFluxConstr
            tmpConstr = optionsMFA.constr.fluxes;
            idTmp = tmpConstr.idRxns;
            for i = 1 : length(idTmp)
                lb.knotFluxesAll(idTmp(i),tmpConstr.idFullKnot(i)) = tmpConstr.lb(i);
                ub.knotFluxesAll (idTmp(i),tmpConstr.idFullKnot(i)) = tmpConstr.ub(i);
            end
        end
        lb.knotFluxesAll = lb.knotFluxesAll/expandRatio;
        ub.knotFluxesAll = ub.knotFluxesAll*expandRatio;
end




end
function idParam = prepIdParam170913(model, expData, optionsMFA, optimType)

field2var(optionsMFA.varSet)
if ~optionsMFA.isUseQPInMH  % 180830修正
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
        lb = optionsMFA.lbMH;
        ub = optionsMFA.ubMH;
    case {'init', 'local'}
        nOptimMetsInitConcs = nNonPoolMets;
        idOptimMetsInitConcs =idNonPoolMets;
        nOptimMetsConcRates = nNonPoolMets;
        idOptimMetsConcRates =idNonPoolMets;
        lb = optionsMFA.lbLocal;
        ub = optionsMFA.ubLocal;
end


%% パラメータのID
nParam = 0;

idParam.initConcs = nParam+(1:nOptimMetsInitConcs);
nParam = nParam + length(idParam.initConcs);

% if isForOptim
%     tmpNKnotsMet = optionsMFA.nKnotsMets(idOptimMets);
%     idParam.knotConcRates = nParam+(1:sum(tmpNKnotsMet+2));
% else
    idParam.knotConcRates = nParam+(1:nOptimMetsConcRates*(nKnots+2));
% end
nParam = nParam + length(idParam.knotConcRates);

% if isForOptim
%     tmpNKnotsRxn = optionsMFA.nKnotsRxns(idIndFluxes);
%     idParam.knotFluxes= nParam+(1:sum(tmpNKnotsRxn+2));
% else
    idParam.knotFluxes = nParam+(1:nIndFluxes * (nKnots+2));
% end
nParam = nParam + length(idParam.knotFluxes);

idParam.knots = nParam + (1:nKnots);
nParam = nParam + length(idParam.knots);

if ~isempty(idComptCorrParam)
    idParam.coefCorrCompt = nParam+(1:nComptCorrParam);
    nParam = nParam + length(idParam.coefCorrCompt);
else
    idParam.coefCorrCompt = [];
end

idParam.nParam =nParam;

%% どのパラメータが最適化されるか
%独立変数の同定のところでやっている。

% isOptParam.initConcs = true(nOptimMetsInitConcs,1);
% loc = ub.initConcs-lb.initConcs<=ub.initConcs*0.01;
% isOptParam.initConcs(loc) = false;
% 
% isOptParam.knotConcRates = true(nOptimMetsConcRates,nKnots+2);
% loc = ub.knotConcRates-lb.knotConcRates<=abs(ub.knotConcRates)*0.01...
%     | ub.knotConcRates-lb.knotConcRates <= 10^-4;
% loc(ismember(idOptimMetsConcRates, model.idRedandantMet),:) = true;% 冗長な代謝物 (NADHとNADの片方など)について
% isOptParam.knotConcRates(loc) = false;
% 
% isOptParam.knotFluxes = true(nIndFluxes,nKnots+2);
% loc = ub.knotFluxes-lb.knotFluxes<=abs(ub.knotFluxes)*0.01;
% isOptParam.knotFluxes(loc) = false;
% 
% isOptParam.knots= true(nKnots,1);
% loc = ub.knots-lb.knots<=abs(ub.knots)*0.01;
% isOptParam.knots(loc) = false;
% 
% isOptParam.coefCorrCompt= true(nComptCorrParam,1);
% loc = ub.coefCorrCompt-lb.coefCorrCompt<=abs(ub.coefCorrCompt)*0.01;
% isOptParam.coefCorrCompt(loc) = false;
% 
% 
% fxnParam2vector = str2func(optionsMFA.fxnName.param2vector);
% isLogTransform = false;
% idParam.isOptParam = fxnParam2vector(...
%     isOptParam, optionsMFA, optimType, isLogTransform); % log10変換しないためここでは必ずlocal;


end

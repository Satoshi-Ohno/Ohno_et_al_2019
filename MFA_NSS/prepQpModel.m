function [qpModel, convertMat, qpModelConstr] = prepQpModel180104(model, expData, optionsMFA,  ...
    optimType, paramInd,isInputIndVars, isInputFullIndVarsConstr)

import packFxnPrepQpModel.*

optimTypeInit = 'init';
convertMat = optionsMFA.convertMat;
if isInputFullIndVarsConstr
    switch optimType
        case {'init', 'metaheuristic'}
                paramMHLin(convertMat.isLogParamMH) = 10.^paramInd(convertMat.isLogParamMH);
        case {'local', 'localNonLinConstr'}
            paramMHLin = [];
        case {'QP'}
            isLogParamMHInd = convertMat.isLogParamMH(optionsMFA.isIndParams.MH);
            paramMHLin(isLogParamMHInd) = 10.^paramInd(isLogParamMHInd);
    end
else
    switch optimType
        case {'init', 'metaheuristic', 'local', 'localNonLinConstr'}
            paramMHLin(convertMat.isLogParamMH) = 10.^paramInd(convertMat.isLogParamMH);
        case {'QP'}
            isLogParamMHInd = convertMat.isLogParamMH(optionsMFA.isIndParams.MH);
            paramMHLin(isLogParamMHInd) = 10.^paramInd(isLogParamMHInd);
    end
end
paramMHLin = colVec(paramMHLin);

%% knots and coefCorrCompt
if isInputIndVars
    switch optimType
        case {'local', 'localNonLinConstr'}
            idParam = optionsMFA.idParamLocal;
            
            % knots
            isInParams = false(1,idParam.nParam);
            isInParams(idParam.knots) = true;
            isInParams = isInParams(optionsMFA.isIndParams.local);
            knots = paramInd(isInParams);
            
            % coefCorrCompt
            isInParams = false(1,idParam.nParam);
            isInParams(idParam.coefCorrCompt) = true;
            isInParams = isInParams(optionsMFA.isIndParams.local);
            coefCorrCompt = paramInd(isInParams);
        case {'metaheuristic', 'init', 'QP'}
            idParam = optionsMFA.idParamMH;
            
            % knots
            isInParams = false(1,idParam.nParam);
            isInParams(idParam.knots) = true;
            isInParams = isInParams(optionsMFA.isIndParams.MH);
            knots = paramMHLin(isInParams);
            
            % coefCorrCompt
            isInParams = false(1,idParam.nParam);
            isInParams(idParam.coefCorrCompt) = true;
            isInParams = isInParams(optionsMFA.isIndParams.MH);
            coefCorrCompt = paramMHLin(isInParams);
%         otherwise
%             keyboard;
    end
else
    idParam = optionsMFA.idParamMH;
    knots = paramMHLin(idParam.knots);
    coefCorrCompt = paramMHLin(idParam.coefCorrCompt);

end
if ~optionsMFA.isAccountMediaDrop && isempty(coefCorrCompt)
    coefCorrCompt = zeros(length(optionsMFA.idParamLocal.coefCorrCompt),1) + ...
        (optionsMFA.lb.coefCorrCompt+optionsMFA.ub.coefCorrCompt)/2;
end

fullKnots = [0, knots', expData.time(end)];

%% Make QP model
fxnPrepConvertMatDepKnots = str2func(optionsMFA.fxnName.prepConvertMatDepKnots);
convertMat = fxnPrepConvertMatDepKnots(model, expData, optionsMFA, fullKnots);
optionsMFA.convertMat = convertMat;

qpModel =  makeQPModelOri(model, expData, optionsMFA, knots, coefCorrCompt, optimTypeInit);

%% modify QP model for constraints in local optimization
if isInputFullIndVarsConstr
    switch optimType
        case {'local'}
            nKnots = optionsMFA.varSet.nKnots;
            
            if nKnots >=1
                idConstrFmincon = [qpModel.idQPConstr.flux1, qpModel.idQPConstr.flux2];
            else 
                idConstrFmincon = [qpModel.idQPConstr.concs, qpModel.idQPConstr.flux1, qpModel.idQPConstr.flux2];
            end
            qpModel.A = qpModel.A(idConstrFmincon, :);
            qpModel.b = qpModel.b(idConstrFmincon);
            
            if nKnots>=1
                AKnotTime = zeros(nKnots+1,idParam.nParam);
                bKnotTime = zeros(nKnots+1,1);
                for k = 1 : nKnots+1
                    switch k
                        case 1
                            AKnotTime(k, idParam.knots(1)) = -1;
                            bKnotTime(k) = -fullKnots(1) - optionsMFA.minKnotTimeDiff;
                        case nKnots+1
                            AKnotTime(k, idParam.knots(nKnots)) = 1;
                            bKnotTime(k) = fullKnots(end) - optionsMFA.minKnotTimeDiff;
                        otherwise
                            AKnotTime(k, idParam.knots(k)) = -1;
                            AKnotTime(k, idParam.knots(k-1)) = -1;
                            bKnotTime(k) =  - optionsMFA.minKnotTimeDiff;
                    end
                end
                qpModel.A = [qpModel.A;sparse(AKnotTime)];
                qpModel.b = [qpModel.b;bKnotTime];
            end
        case {'localNonLinConstr'}
            idConstrFmincon = [qpModel.idQPConstr.concs];
            qpModel.A = qpModel.A(idConstrFmincon, :);
            qpModel.b = qpModel.b(idConstrFmincon);
    end
end

%% Add mass balance constraints
qpModel = modifyQPModelMassBalance(model, expData, optionsMFA, qpModel, optimType);

%% Add constrainsf for reaction dependent number of time intervals
if optionsMFA.isRxnSpecificNKnots
    qpModel = modifyQPModelConstrRxnSpecificNKnots(model, expData, optionsMFA, ...
        qpModel, optimTypeInit);
end

%% Add constrain for paramMH
switch optimType
    case {'init', 'QP'}
        qpModel = modifyQPModelWithParamMH(model, expData, optionsMFA, ...
            qpModel, optimType, paramInd);
    case {'local', 'metaheuristic','localNonLinConstr'}
end
qpModelConstr = qpModel;

%% Remove lb and ub and add linier constrains to A or Aeq
qpModel = modifyQPModelWithLbUb(model, expData, optionsMFA, qpModel, optimType);

%% Identify independent parameters from Aeq
fxnSelectIndVars = str2func(optionsMFA.fxnName.selectIndVars);
if isInputFullIndVarsConstr
    switch optimType
        case {'local', 'localNonLinConstr'}
            idInputIndVars = find(optionsMFA.isIndParams.local);
            idInputRedundantConstr = find(optionsMFA.isRedundantConstrAeq.local);
        case {'init'}
            idInputIndVars = find(optionsMFA.isIndParams.init);
            idInputRedundantConstr = find(optionsMFA.isRedundantConstrAeq.init);
        otherwise
            idInputIndVars = find(optionsMFA.isIndParams.QP);
            idInputRedundantConstr = find(optionsMFA.isRedundantConstrAeq.QP);
    end
    if isempty(full(qpModel.Aeq))
        idIndVars = (1:size(qpModel.A,2))';
        idRedundantConstr = [];
        Cx = speye(size(qpModel.A,2));
        C0 = [];
    else
        [idIndVars, idRedundantConstr, Cx, C0] = ...
            fxnSelectIndVars(full(qpModel.Aeq), isInputFullIndVarsConstr, idInputIndVars, idInputRedundantConstr);
    end
    
     
else 
    if isempty(full(qpModel.Aeq))
        idIndVars = (1:size(qpModel.A,2))';
        idRedundantConstr = [];
        Cx = speye(size(qpModel.A,2));
        C0 = [];
    else
        [idIndVars, idRedundantConstr, Cx, C0] = ...
            fxnSelectIndVars(full(qpModel.Aeq), isInputFullIndVarsConstr, []);
    end
    
        switch optimType
            case {'metaheuristic'}
                idParamLocal = optionsMFA.idParamLocal;
                nNonPoolMets = optionsMFA.varSet.nNonPoolMets;
                nKnots = optionsMFA.varSet.nKnots;
                nIndFluxes = optionsMFA.varSet.nIndFluxes;
                
                idParamKnotConcRates = reshape(idParamLocal.knotConcRates, nNonPoolMets, nKnots+2);
                idParamKnotFluxes = reshape(idParamLocal.knotFluxes, nIndFluxes, nKnots+2);
                orderVar{1} = [...
                    idParamLocal.initConcs,...
                    idParamKnotConcRates(:,1)',...
                    idParamKnotConcRates(:,end)',...
                    idParamKnotFluxes(:,1)',...
                    idParamKnotFluxes(:,end)',...
                    idParamLocal.knots,...
                    idParamLocal.coefCorrCompt,...
                    ];
                if nKnots >= 1
                    orderVar{1} = [orderVar{1},...
                        reshape(idParamKnotConcRates(:,2:end-1), nNonPoolMets*nKnots, 1)',...
                        reshape(idParamKnotFluxes(:,2:end-1), nIndFluxes*nKnots, 1)',...
                        ];
                end
                if ~isempty(qpModel.Aeq)
                    orderVar{1} =  orderVar{1}(end:-1:1);
                    orderConstr{1} = 1:size(qpModel.Aeq,1);
                    for i = 1 : length(orderVar)
                        [tmpIdIndVars, tmpIdRedundantConstr, Cx, C0] = ...
                            fxnSelectIndVars(full(qpModel.Aeq(orderConstr{i}, orderVar{i})), isInputFullIndVarsConstr, []);
                        isIndVars = false(size(qpModel.Aeq, 2),1);
                        isIndVars(tmpIdIndVars) = true;
                        [~, idSort] = sort(orderVar{i});
                        isIndVars = isIndVars(idSort);
                        idIndVarsTest = find(isIndVars);
                        
                        isRedundantConstr = false(size(qpModel.Aeq, 1),1);
                        isRedundantConstr(tmpIdRedundantConstr) = true;
                        [~, idSort] = sort(orderConstr{i});
                        isRedundantConstr = isRedundantConstr(idSort);
                        idRedundantConstrTest = find(isRedundantConstr);
                        idIndVars = idIndVarsTest;
                        idRedundantConstr = idRedundantConstrTest;
                    end
                end
        end
end
qpModel.idIndVars = idIndVars;
qpModel.idRedundantConstr = idRedundantConstr;

%% Remove Aeq and 
qpModel = modifyQPModelRmAeq(model, expData, optionsMFA, qpModel, Cx, C0);


end
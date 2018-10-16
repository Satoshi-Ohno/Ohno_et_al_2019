%% Make initial parameters
function [params, scores, nExpData] = createInitParams171226(model, expData,optionsMFA)

field2var(optionsMFA.varSet)

popSize = optionsMFA.optionsOptimMH.popSize;
idParamMH = optionsMFA.idParamMH;
nParamMH = idParamMH.nParam;
nParamMHInd = nnz(optionsMFA.isIndParams.MH);
lbInit = optionsMFA.lbInit;
ubInit = optionsMFA.ubInit;
optimType = 'init';

%% Set objective function
fxnSolveMDVDynamics = str2func(optionsMFA.fxnName.solveMDVDynamics);
objfunMH =@(params) fxnSolveMDVDynamics(...
    params, model, expData, optionsMFA, 'metaheuristic');
objfunLocal =@(params) fxnSolveMDVDynamics(...
    params, model, expData, optionsMFA, 'local');


%% Make initial parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Making initial parameters')
p = 1;
initPopSize = popSize;
params = zeros(nParamMHInd,initPopSize);
scores = zeros(1,initPopSize);
nFailureParamConstr = 0;
nFailureODE = 0;
ii  = 1;
while p <= initPopSize
    nFailureLimit = max([10^5, initPopSize*10^4]);
    if nFailureParamConstr >= nFailureLimit;
        error(['Initial parameter sets were not found within ' num2str(nFailureLimit) ' trials'])
    end
    
    tmpParamMH = zeros(nParamMH, 1);

    %% Initial knot time and coefCorrCompt
    % knot time
    log10Range = log10(ubInit.knots)-log10(lbInit.knots);
    while true
        knots = 10.^(log10Range.*rand(nKnots, 1)+log10(lbInit.knots));
        knots = sort(knots);
        fullKnots = [0, knots', expData.time(end)];
        if all(diff(fullKnots) >= optionsMFA.minKnotTimeDiff)
            break
        end
    end

    % coefCorrCompt
    if ~isempty(idComptCorrParam)
        log10Range = log10(ubInit.coefCorrCompt)-log10(lbInit.coefCorrCompt);
        coefCorrCompt = 10.^(log10Range.*rand(nComptCorrParam, 1)+log10(lbInit.coefCorrCompt));
    else
        mediaInfo = optionsMFA.mediaInfo;
        coefCorrCompt = 1/mediaInfo.mediaVol*mediaInfo.cellAmount;
    end
    
    tmpParamMH(idParamMH.knots) = log10(knots);
    if ~isempty(idParamMH.coefCorrCompt)
        tmpParamMH(idParamMH.coefCorrCompt) = log10(coefCorrCompt);
    end
    
    %% Make QP model
    isInputIndVars = false;
    isInputFullIndVarsConstr= true;
    fxnPrepQpModel = str2func(optionsMFA.fxnName.prepQpModel);
    [qpModel, qpModelConstr] = fxnPrepQpModel(model, expData, optionsMFA,  ...
        optimType, tmpParamMH,isInputIndVars, isInputFullIndVarsConstr);

    %% Modify objective function of QP
    isMinVar = false;
    fxnModifyQPModelObjFun = str2func(optionsMFA.fxnName.modifyQPModelObjFun);
    qpModel = fxnModifyQPModelObjFun(model, expData, optionsMFA, qpModel, isMinVar);
    
    %% Solve QP to fit metabolite concentrations
    [qpSol, qpScore,exitFlag] = solveQP(qpModel);
    tmpParamInitInd = qpSol;


    
    %% paramLocal to paramMH
    tmpParamLocalFull = qpModel.L * tmpParamInitInd + qpModel.m;
    tmpParamMHFull = paramLocal2ParamMH(model, expData, optionsMFA, tmpParamLocalFull);
    if optionsMFA.isUseQPInMH
        tmpParamMH = tmpParamMHFull(optionsMFA.isIndParams.MH);
    else
        tmpParamMH = tmpParamMHFull(optionsMFA.isIndParams.local);
    end
        
    %% Check feasibility
    nFailureODE = nFailureODE + 1;    
    [tmpScore, ~, sol] = objfunMH(tmpParamMH);
   
    if tmpScore < 10^12
        params(:,p) = tmpParamMH;
        scores(p) = tmpScore;
        if rem(p,ceil(popSize/5)) == 0
            disp(['Initial param#' num2str(p)])
        end
        p = p + 1;
        nFailureODE = 0;
        nFailureParamConstr = 0;
        continue
    end
    if nFailureODE >= 10^4
        nFailureODE = 0;
        break
    end
end

%% Sort paramter sets according to the scores
[scores, idSort] = sort(scores);
scores=scores(1:popSize);
params= params(:,idSort(1:popSize));
nExpData = sol.nExpData;

end


%% Make initial parameters
function [params, scores, nExpData] = createInitParams171226(model, expData,optionsMFA)

field2var(optionsMFA.varSet)

initPopSize = optionsMFA.optionsOptimMH.initPopSize;
idParamMH = optionsMFA.idParamMH;
nParamMH = idParamMH.nParam;
nParamMHInd = nnz(optionsMFA.isIndParams.MH);
lbInit = optionsMFA.lbInit;
ubInit = optionsMFA.ubInit;
optimType = 'init';

%% Set objective function
objfunMH =@(params) solveMDVDynamics(...
    params, model, expData, optionsMFA, 'metaheuristic');
objfunLocal =@(params) solveMDVDynamics(...
    params, model, expData, optionsMFA, 'local');


%% Make initial parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Making initial parameters')
p = 1;
params = zeros(nParamMHInd,initPopSize);
scores = zeros(1,initPopSize);
nFailureParamConstr = 0;
nFailureODE = 0;
ii  = 1;
while p <= initPopSize
    nFailureLimit = max([10^5, initPopSize*10^4]);
    if nFailureParamConstr >= nFailureLimit
        error(['Initial parameter sets were not found within ' num2str(nFailureLimit) ' trials'])
    end
    
    tmpParamMH = zeros(nParamMH, 1);

    %% Initial switch time 
    % switch time
    log10Range = log10(ubInit.switchTimes)-log10(lbInit.switchTimes);
    while true
        switchTimes = 10.^(log10Range.*rand(nSwitchTimes, 1)+log10(lbInit.switchTimes));
        switchTimes = sort(switchTimes);
        fullSwitchTimes = [0, switchTimes', expData.time(end)];
        if all(diff(fullSwitchTimes) >= optionsMFA.minKnotTimeDiff)
            break
        end
    end
    
    tmpParamMH(idParamMH.switchTimes) = log10(switchTimes);
    
    %% Metabolite concentration
    if optionsMFA.isUseConcAsParam
        log10Range = log10(ubInit.switchTimeConcs)-log10(lbInit.switchTimeConcs);
        concs = 10.^(log10Range.*rand(size(ubInit.switchTimeConcs))+log10(lbInit.switchTimeConcs));
        if optionsMFA.isSameInitInsCtrl
            concs(nNonPoolMets/2 + (1:nNonPoolMets/2),1) = concs(1:nNonPoolMets/2,1);
        end
        tmpParamMH(idParamMH.concs) = log10(concs(:));
    end
    
    %% Make QP model
    isInputIndVars = false;
    isInputFullIndVarsConstr= true;
    [qpModel, qpModelConstr] = prepQpModel(model, expData, optionsMFA,  ...
        optimType, tmpParamMH,isInputIndVars, isInputFullIndVarsConstr);

    %% Modify objective function of QP
    isMinVar = false;
    qpModel = modifyQPModelObjFun(model, expData, optionsMFA, qpModel, isMinVar);
    
    %% Solve QP to fit metabolite concentrations
    [qpSol, qpScore,exitFlag] = solveQP(qpModel);
    tmpParamInitInd = qpSol;
    if isempty(qpSol)
        continue
    end
    if exitFlag == -2
        continue        
    end
    
    %% Sample parameters uniformly under linear constraints
    try
        tmpParamInitInd = smplParamsUnderLinConstr(model, expData, optionsMFA, qpModel, tmpParamInitInd);
    catch err
        tmpParamInitInd = [];
    end
    if isempty(tmpParamInitInd)
        nFailureParamConstr = nFailureParamConstr + 1;
        continue
    end
    
    %% paramLocal to paramMH
    tmpParamLocalFull = qpModel.L * tmpParamInitInd + qpModel.m;
    optionsMFA.transformMat = prepTransformMat(model, expData, optionsMFA, fullSwitchTimes);
    tmpParamStruct = vector2paramStruct(tmpParamLocalFull, model, expData, optionsMFA);
    testConcsExpTime = optionsMFA.transformMat.param2ConcExpTime * tmpParamLocalFull;
    testConcsExpTime = reshape(testConcsExpTime, nExpTime, optionsMFA.varSet.nEvalConcMets)';
    testConcsQP = optionsMFA.transformMat.param2ConcQP* tmpParamLocalFull;
    testConcsQP = reshape(testConcsQP, optionsMFA.nStepTimeConstrQP,nNonPoolMets)';
    
    tmpParamMHFull = paramLocal2ParamMH(model, expData, optionsMFA, tmpParamLocalFull);
    if optionsMFA.isUseQPInMH
        tmpParamMH = tmpParamMHFull(optionsMFA.isIndParams.MH);
    else
        tmpParamMH = tmpParamMHFull(optionsMFA.isIndParams.local);
    end
    if ~ isreal(tmpParamMH)
        continue
    end
        
    %% Check feasibility
    nFailureODE = nFailureODE + 1;    
    
    [tmpScore, ~, sol] = objfunMH(tmpParamMH);
   
    if tmpScore < 10^12
        params(:,p) = tmpParamMH;
        scores(p) = tmpScore;
        if rem(p,ceil(initPopSize/5)) == 0
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
params= params(:,idSort);
nExpData = sol.nExpData;

end


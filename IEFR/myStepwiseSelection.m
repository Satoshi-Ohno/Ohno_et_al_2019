function     [solFinal, solInit] = myStepwiseFit(Xall, objfun, options)
%% stepwise model selection

[nData, nParam] = size(Xall);

isConstrParam = options.ub==options.lb;

isUseVar = false(1, nParam);
isUseVar(isConstrParam) = options.ub(isConstrParam);
isUseVar(1) = true;

nStep = 0;
flgStop = false;
isAddParam = true;

%% step i: Set initial model
isUseVarMat = isUseVar;
nChangeUseVar = 1;

% evaluate
score = zeros(nChangeUseVar,1);
for i =  1: nChangeUseVar
    isUseVar =isUseVarMat(i,:);
    [score(i), sol(i)] = objfun(isUseVar);
end
[scoreMin, idMin] = min(score);
solMin = sol(idMin);
isUseVar = isUseVarMat(idMin,:);

scoreBest = scoreMin;
solBest = solMin;
isUseVarBest = isUseVar;

% save current step
nStep = nStep + 1;
[~,idSort] = sort(score);
tmpHistory.nStep = nStep;
tmpHistory.isAddParam = isAddParam;
tmpHistory.scoreMin = scoreMin;
tmpHistory.isUseVarMin = isUseVar;
tmpHistory.score = score(idSort);
tmpHistory.isUseVarMat = isUseVarMat(idSort,:);
tmpHistory.updated = true;
history(nStep) = tmpHistory;

% display current step
disp(' ')
disp(['# of step: ' num2str(nStep) ...
    ',   # of param: '  num2str(nnz(isUseVar)) ...
    ',   score: ' num2str(scoreMin)])

% save initial model
solInit.score = scoreBest;
solInit.sol = solBest;
solInit.isUseVar = isUseVarBest;
solInit.history = history;

while true
    %% step ii: Evaluate models added with each of parameters
    isAddParam = true;
    idChangeUseVar = find(isUseVar==0 & options.ub'==1);
    if isempty(idChangeUseVar)
        flgStop = true;
    else
        nChangeUseVar = length(idChangeUseVar);
        isUseVarMat = repmat(isUseVar, nChangeUseVar,1);
        isUseVarMat(:,idChangeUseVar) = isUseVarMat(:,idChangeUseVar) + eye(nChangeUseVar);
    end
    
    if ~flgStop
        % evaluate
        score = zeros(nChangeUseVar,1);
        for i =  1: nChangeUseVar
            isUseVar =isUseVarMat(i,:);
            [score(i), sol(i)] = objfun(isUseVar);
        end
        [scoreMin, idMin] = min(score);
        solMin = sol(idMin);
        isUseVar = isUseVarMat(idMin,:);
        
    end
    
    %% step iii: Add one parameter
    if ~flgStop
        % save current step
        nStep = nStep + 1;
        [~,idSort] = sort(score);
        tmpHistory.nStep = nStep;
        tmpHistory.isAddParam = isAddParam;
        tmpHistory.scoreMin = scoreMin;
        tmpHistory.isUseVarMin = isUseVar;
        tmpHistory.score = score(idSort);
        tmpHistory.isUseVarMat = isUseVarMat(idSort,:);
        
        % display current step
        disp(' ')
        disp(['# of step: ' num2str(nStep) ...
            ',   # of param: '  num2str(nnz(isUseVar)) ...
            ',   score: ' num2str(scoreMin)])
        
        % Update or reject
        if scoreMin < scoreBest
            disp('Updated.')
            tmpHistory.updated = true;
            scoreBest = scoreMin;
            solBest = solMin;
            isUseVarBest = isUseVar;
        else
            disp('Rejected.')
            tmpHistory.updated = false;
            isUseVar=isUseVarBest;
            flgStop = true;
        end
        history(nStep) = tmpHistory;
    end
    
    %% step iv: Evaluate models without each of parameters
    if ~flgStop
        isAddParam = false;
        flgGoToStep2 = false;
        while true
            idChangeUseVar = find(isUseVar==1 & options.lb'==0);
            if isempty(idChangeUseVar)
                flgGoToStep2 = true;
            else
                nChangeUseVar = length(idChangeUseVar);
                isUseVarMat = repmat(isUseVar, nChangeUseVar,1);
                isUseVarMat(:,idChangeUseVar) = isUseVarMat(:,idChangeUseVar) - eye(nChangeUseVar);
            end
            
            if ~flgGoToStep2
                % evaluate
                score = zeros(nChangeUseVar,1);
                for i =  1: nChangeUseVar
                    isUseVar =isUseVarMat(i,:);
                    [score(i), sol(i)] = objfun(isUseVar);
                end
                [scoreMin, idMin] = min(score);
                solMin = sol(idMin);
                isUseVar = isUseVarMat(idMin,:);
                
            end
            
            %% step v: Remove one parameter
            if ~flgGoToStep2
                % save current step
                nStep = nStep + 1;
                [~,idSort] = sort(score);
                tmpHistory.nStep = nStep;
                tmpHistory.isAddParam = isAddParam;
                tmpHistory.scoreMin = scoreMin;
                tmpHistory.isUseVarMin = isUseVar;
                tmpHistory.score = score(idSort);
                tmpHistory.isUseVarMat = isUseVarMat(idSort,:);
                
                % display current step
                disp(' ')
                disp(['# of step: ' num2str(nStep) ...
                    ',   # of param: '  num2str(nnz(isUseVar)) ...
                    ',   score: ' num2str(scoreMin)])
                
                % update or reject
                if scoreMin < scoreBest
                    disp('Updated.')
                    tmpHistory.updated = true;
                    scoreBest = scoreMin;
                    solBest = solMin;
                    isUseVarBest = isUseVar;
                else
                    disp('Rejected.')
                    tmpHistory.updated = false;
                    isUseVar=isUseVarBest;
                    flgGoToStep2 = true;
                end
                history(nStep) = tmpHistory;
            end
            if flgGoToStep2
                break
            end
        end
    end
    
    %% step vi: Stop
    if flgStop
        disp(' ')
        disp('Finished.')
        disp(['# of step: ' num2str(nStep) ...
            ',   # of param: '  num2str(nnz(isUseVarBest)) ...
            ',   score: ' num2str(scoreBest)])
        break
    end
    
end

%% final model
solFinal.score = scoreBest;
solFinal.sol = solBest;
solFinal.isUseVar = isUseVarBest;
solFinal.history = history;
end


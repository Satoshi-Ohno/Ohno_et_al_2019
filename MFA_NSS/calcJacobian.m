%% Calculate Jacobian of objective function
function      [jacobian1, jacobian2, paramDiffTypes] = calcJacobian(objfunRes, param, lb, ub, optionsCalcJacobian)

field2var(optionsCalcJacobian)

nParam = length(param);
switch evalFxnType
    case {'obj'}
        resOri1 = objfunRes(param);
        resOri2 = [];
    case {'nonlcon'}
        [resOri1, resOri2]= objfunRes(param);
end
optionsCalcJacobian.resOri1 = resOri1;
optionsCalcJacobian.resOri2 = resOri2;

jacobian1 = [];
jacobian2 = [];

%% jacobian‚ÌŒvŽZ
nData1 = length(resOri1);
diffResMat1 = zeros(nData1, nParam);
diffParamMat1 = zeros(nData1, nParam);
switch evalFxnType
    case {'obj'}
    case {'nonlcon'}
        if isempty(resOri2)
        else
            nData2 = length(resOri2);
            diffResMat2 = zeros(nData2, nParam);
            diffParamMat2 = zeros(nData2, nParam);
        end
end


diffTypeAll = optionsCalcJacobian.diffType;
paramDiffTypes = zeros(nParam,1);
for p = 1 : nParam
    isSuccess = false;
        for i = 1 : length(diffTypeAll)
            optionsCalcJacobian.diffType = diffTypeAll{i};
            try
                [tmpDiffRes1, tmpDiffParam1, tmpDiffRes2, tmpDiffParam2] = ...
                    calcDiffRes(objfunRes, param, lb, ub, optionsCalcJacobian, p);
                isSuccess = true;
            catch err
            end
            if isSuccess
                paramDiffTypes(p) = i;
                break
            end
        end
        if ~isSuccess
            paramDiffTypes(p) = 0;
            tmpDiffRes1 = inf(nData1, 1);
            tmpDiffParam1 = zeros(nData1, 1)+0.1;
            if strcmp({evalFxnType}, {'nonlcon'})
                if isempty(resOri2)
                else
                    tmpDiffRes2 = inf(nData2, 1);
                    tmpDiffParam2 = zeros(nData2, 1)+0.1;
                end
            end
%             return
        end
        
    diffResMat1(:,p) = tmpDiffRes1;
    diffParamMat1(:,p) = tmpDiffParam1;
    if strcmp({evalFxnType}, {'nonlcon'}) 
        if isempty(resOri2)
        else
            diffResMat2(:,p) = tmpDiffRes2;
            diffParamMat2(:,p) = tmpDiffParam2;
        end
    end
end

jacobian1 = diffResMat1 ./ diffParamMat1;
if strcmp({evalFxnType}, {'nonlcon'}) 
    if isempty(resOri2)
        jacobian2 = [];
    else
        jacobian2 = diffResMat2 ./ diffParamMat2;
    end
end

end

%%  Calculate difference
function [tmpDiffRes1, tmpDiffParam1, tmpDiffRes2, tmpDiffParam2] = ...
    calcDiffRes(objfunRes, param, lb, ub, optionsCalcJacobian, p)

diffType = optionsCalcJacobian.diffType;
evalFxnType = optionsCalcJacobian.evalFxnType;
initDeltaRatio = optionsCalcJacobian.initDeltaRatio;
resOri1 = optionsCalcJacobian.resOri1;
resOri2 = optionsCalcJacobian.resOri2;

nParam = length(param);
nData1 = length(resOri1);
nData2 = length(resOri2);

tmpDiffRes1=[]; tmpDiffParam1=[];  tmpDiffRes2=[]; tmpDiffParam2=[];

deltaRatio = initDeltaRatio;


nWhileLoop = 0;
tmpSign = -1;
minAbsBound = min(abs([lb(p), ub(p)]));
while true
    nWhileLoop = nWhileLoop + 1;

    tmpParam = [param, param];
    tmpRes1 = [resOri1, resOri1];
    tmpRes2 = [resOri2, resOri2];
    switch diffType
        case {'forward'}
            i = 2;
            tmpParam(p,i) = tmpParam(p,i) + tmpSign * abs(tmpParam(p,i))*deltaRatio;
            if tmpParam(p,i) > lb(p)  && tmpParam(p,i) < ub(p)
                switch evalFxnType
                    case {'obj'}
                        tmpRes1(:,i) = objfunRes(tmpParam(:,i));
                    case {'nonlcon'}
                        [tmpRes1(:,i), tmpRes2(:,i)] = objfunRes(tmpParam(:,i));
                end
            end
        case {'central'}
            for i = 1 : 2
                if i== 1
                    tmpSign = -1;
                else
                    tmpSign = +1;
                end
                tmpParam(p,i) = tmpParam(p,i) + tmpSign * abs(tmpParam(p,i))*deltaRatio;
                if tmpParam(p,i) > lb(p)  && tmpParam(p,i) < ub(p)
                    switch evalFxnType
                        case {'obj'}
                            tmpRes1(:,i) = objfunRes(tmpParam(:,i));
                        case {'nonlcon'}
                            [tmpRes1(:,i), tmpRes2(:,i)] = objfunRes(tmpParam(:,i));
                    end
                end
            end
    end
    
    % Stop flag
    flgStop = true;
    if any(tmpParam(p,:) < lb(p))
        flgStop = false;
    end
    if any(tmpParam(p,:) > ub(p))
        flgStop = false;
    end
    if strcmp({evalFxnType}, {'obj'})
        if any(sum(tmpRes1.^2,1) >= 10^9)
            flgStop = false;
        end
    end
    
    if nWhileLoop >= 30
        error(['Parameter #' num2str(p) 'cannnot be perturbed at ' num2str(deltaRatio*100) '%'])
    end
    if flgStop
        break
    end
    switch diffType
        case {'cetnral'}
            deltaRatio = 0.5*deltaRatio;
        case {'forward'}
            if tmpSign==-1
                tmpSign = - tmpSign;
            else
                tmpSign = - tmpSign;
                deltaRatio = 0.5*deltaRatio;
            end
    end
end  % end while

tmpDiffRes1 = diff(tmpRes1, [], 2);
tmpDiffParam1 = zeros(nData1, 1) + diff(tmpParam(p,:));
switch evalFxnType
    case {'obj'}
        tmpDiffRes2 = [];
        tmpDiffParam2 = [];
    case {'nonlcon'}
        tmpDiffRes1 = diff(tmpRes2, [], 2);
        tmpDiffParam1 = zeros(nData2, 1) + diff(tmpParam(p,:));
end

end



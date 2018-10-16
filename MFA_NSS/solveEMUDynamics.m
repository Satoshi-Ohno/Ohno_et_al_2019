%% Simulation of EMU
function output =  solveEMUDynamics170809(paramStruct, model, expData, optionsMFA, preSol)

field2var(optionsMFA.varSet)

querryTime = paramStruct.querryTime;
flagFailure = 0;

nOdeStep = [];
for s = 1 : length(model.emuNetwork)
    emuNetwork = model.emuNetwork(s);
    arrayMakeA = emuNetwork.arrayMakeA;
    arrayMakeB = emuNetwork.arrayMakeB;  
    emuS = emuNetwork.emuS;
    if isempty(arrayMakeA)
        continue
    end
    if isempty(arrayMakeA.array{1})
        continue
    end
    unknownEmu = emuNetwork.arrayMakeA.emu;
    knownEmu = emuNetwork.arrayMakeB.emu;
    nUnkEmu = length(unknownEmu);
    nKnownEmu = length(knownEmu);    
    emu = [unknownEmu;knownEmu];
    clear emuNetwork
    
    %%  Calculate A, B, and S
    A0 = zeros(nUnkEmu);
    A1 = zeros(nUnkEmu);
    for i = 1 : nRxns
        A0 = A0 + arrayMakeA.array{i} * paramStruct.pCoefFluxes(i,1);
        A1 = A1 + arrayMakeA.array{i} * paramStruct.pCoefFluxes(i,2);
    end
    B0 = zeros(nUnkEmu, nKnownEmu);
    B1 = zeros(nUnkEmu, nKnownEmu);
    for i = 1 : nRxns
        B0 = B0 + arrayMakeB.array{i} * paramStruct.pCoefFluxes(i,1);
        B1 = B1 + arrayMakeB.array{i} * paramStruct.pCoefFluxes(i,2);
    end
    Sv0 = emuS*paramStruct.pCoefFluxes(:,1);
    Sv1 = emuS*paramStruct.pCoefFluxes(:,2);
    
    paramStruct.A0 = A0;
    paramStruct.A1 = A1;
    paramStruct.B0 = B0;
    paramStruct.B1 = B1;
    paramStruct.Sv0 = Sv0;
    paramStruct.Sv1 = Sv1;
%     clear A B Sv
    
    %% polinomial coefficient of metabolite concentrations
    paramStruct.tmpPCoefConcs = paramStruct.pCoefConcs([arrayMakeA.idMet],:);

    %% Known EMU
    
    result(s).emu = emu;
    result(s).unknownEmu = unknownEmu;
    result(s).unknownEmuDynamics = [];
    result(s).unknownEmuTime = [];
    result(s).knownEmu = knownEmu;
    result(s).knownEmuDynamics = [];
    result(s).knownEmuTime = [];
    
    knownEmuDynamics = cell(nKnownEmu, s+1);
    for e = 1 : nKnownEmu
        if arrayMakeB.isSingleEmu(e)
            poolEmuDynamics = makePoolEmuDynamics(model, querryTime, emu(nUnkEmu+e), s, expData,optionsMFA.isCorrectNatAbundanceSim);
            knownEmuDynamics(e,:) = poolEmuDynamics;
        else   
            refEmuSize = arrayMakeB.sizeEmu{e};
            tmpEmuDynamics = cell(2, 10);
            for i = 1 : 2
                idRefEmu = find(strcmp(emu{nUnkEmu+e}(i), result(refEmuSize(i)).unknownEmu)); 
                if isempty(idRefEmu)
                    poolEmuDynamics = makePoolEmuDynamics(model, querryTime, emu{nUnkEmu+e}(i), refEmuSize(i), expData,optionsMFA.isCorrectNatAbundanceSim);
                    tmpEmuDynamics(i,1:refEmuSize(i)+1) = poolEmuDynamics;
                else
                    for k = 1 : refEmuSize(i)+1
                        refEmuDynamics = result(refEmuSize(i)).unknownEmuDynamics{k}(idRefEmu,:);
                        tmpEmuDynamics{i,k} = refEmuDynamics;
                    end
                end
            end
            massMatrix = zeros(refEmuSize(1)+1, refEmuSize(2)+1);
            massMatrix = massMatrix + repmat([0:refEmuSize(1)]', 1, refEmuSize(2)+1);
            massMatrix = massMatrix + repmat([0:refEmuSize(2)] , refEmuSize(1)+1, 1);
            for k = 1 : s+1
                [idr,idc] = find(massMatrix == k-1);
                knownEmuDynamics{e,k} = zeros(1,length(querryTime));
                for i = 1 : length(idr)
                knownEmuDynamics{e,k} = knownEmuDynamics{e,k} + ...
                    tmpEmuDynamics{1, idr(i)} .*  tmpEmuDynamics{2, idc(i)};
                end
            end            
        end
        result(s).knownEmuDynamics = knownEmuDynamics;
    end
    result(s).knownEmuTime = querryTime;
    
    paramStruct.Y = result(s).knownEmuDynamics;
    paramStruct.t = result(s).knownEmuTime;
    clear knownEmuDynamics
        
    %% Simulate unknown EMU
    result(s).unknownEmuDynamics = cell(1,s+1);
    result(s).unknownEmuTime = cell(1,s+1);
    
    tmpOutput = simulateEmuDynamics(...
        model, querryTime, paramStruct, s, expData.nat13CRatio, preSol, optionsMFA);

   sumEmuOverM0 = zeros(nUnkEmu, length(tmpOutput.time));
    for k = 2 : s+1
        loc = (k-2)*nUnkEmu+(1:nUnkEmu);
        result(s).unknownEmuDynamics{k} =tmpOutput.emu(loc,:);
        sumEmuOverM0 = sumEmuOverM0 + tmpOutput.emu(loc,:);
    end
        
    if length(tmpOutput.time)~=length(querryTime)
        if s == 1
            querryTime=tmpOutput.time;
        else
            disp('   Time for ODE is different between input and output.')
        end
    end
    result(s).unknownEmuTime =querryTime;
        
    maxEmu = 1;
    minEmu = 0;
    idMax = 0;
    for k = 2 : s+1
        tmpEmuDynamics = result(s).unknownEmuDynamics{k};
        maxEmu = max([maxEmu, reshape(tmpEmuDynamics, 1, numel(tmpEmuDynamics))]);
        if maxEmu > 1+optionsMFA.tolEmuError
            [maxVal, idMax] = mmax(tmpEmuDynamics);
            break
        end
    end
    idMin = 0;
    for k = 2 : s+1
        tmpEmuDynamics = result(s).unknownEmuDynamics{k};
        minEmu = min([minEmu, reshape(tmpEmuDynamics, 1, numel(tmpEmuDynamics))]);
        if minEmu < 0-optionsMFA.tolEmuError
            [minVal, idMin] = mmin(tmpEmuDynamics);
            break
        end
    end
    if maxEmu > 1+optionsMFA.tolEmuError || minEmu < 0-optionsMFA.tolEmuError
        disp(['      ODE is not stable (' unknownEmu{max([idMax(1), idMin(1)])} '). Max emu=' num2str(maxEmu,3), ', Min emu=' num2str(minEmu,3) ', EMU size=' num2str(s)])
        flagFailure = 1;
        break
    elseif maxEmu > 1 || minEmu < 0
    end
    
    for k = 2 : s+1
        tmpEmuDynamics = result(s).unknownEmuDynamics{k};        
        tmpEmuDynamics(tmpEmuDynamics <= 0) = 0;
        tmpEmuDynamics(tmpEmuDynamics >= 1) = 1;
        result(s).unknownEmuDynamics{k} = tmpEmuDynamics ;        
    end
    result(s).unknownEmuDynamics{1} =1 - sumEmuOverM0;
        
end

%% Orgranise results
if flagFailure == 1
    output.timeSim = querryTime;
    output.MDVsSim = [];
    output.EMUresult = [];
    output.concs = [];
    output.knotFluxes = paramStruct.knotFluxes;
    output.nOdeStep = nOdeStep;
    return
end
MDVsSim = cell(nMets,1);
for i = 1 : nMets
    sizeEmu = model.metInfo.nCarbonMets(i);
    tmpMDVsSim = zeros(sizeEmu+1, length(querryTime));
    
    tmpEMUNumber = 0;
    for k = 1 : sizeEmu
        tmpEMUNumber = tmpEMUNumber *10+k;
    end
    emu = [model.mets{i} num2str(tmpEMUNumber )];
    idEmu = find(strcmp(result(sizeEmu).unknownEmu, {emu}));
    if isempty(idEmu)
        switch model.metInfo.massType(i)
            case 0
            case 1
                tmpMDVsSim(end,:) = 1;
            otherwise
                tmpMDVsSim(1,:) = 1;
        end
    else
        for j = 1 : sizeEmu+1
            refEmuDynamics = result(sizeEmu).unknownEmuDynamics{j}(idEmu,:);
            tmpMDVsSim(j,:) =refEmuDynamics;
        end
    end
    MDVsSim{i} = tmpMDVsSim;
end

for s1 = 1 : length(result)
    for s2 = 1 : length(result(s1).unknownEmuDynamics)
        result(s1).unknownEmuDynamics{s2} = result(s1).unknownEmuDynamics{s2}(:,[1,end]);
    end
end
result = rmfield(result, {'knownEmuDynamics', 'knownEmuTime', 'unknownEmuTime'});

output.timeSim = querryTime;
output.MDVsSim = MDVsSim;
output.EMUresult = result;
output.concs = paramStruct.pCoefConcs * [querryTime.^0;querryTime.^1;querryTime.^2];
output.knotFluxes = paramStruct.knotFluxes;
output.nOdeStep = nOdeStep;

end

%% Simulate unknown EMU
function output = simulateEmuDynamics(...
    model, querryTime, params, sizeEmu, nat13CRatio, preSol, optionsMFA)
% Note that M+0 mass isotopomer are not simulated
 
optionsODE = optionsMFA.optionsODE;

nUnkEmu = size(params.A0,2);
nKnownEmu = size(params.B0,2);
initX = zeros(nUnkEmu,sizeEmu) ;
if isempty(preSol)
    for k = 2 : sizeEmu+1
        initX(:,k-1) = zeros(nUnkEmu,1) + ...
            nchoosekScalar(sizeEmu,k-1) * nat13CRatio(2)^(k-1) * nat13CRatio(1)^(sizeEmu-k+1);
    end
else
    for k = 2 : sizeEmu+1
        initX(:,k-1) = [preSol.EMUresult(sizeEmu).unknownEmuDynamics{k}(:,end)];
    end
end
initX = reshape(initX, nUnkEmu*(sizeEmu),1);

odeFunc = str2func(optionsODE.solver);

paramsODE.A0 = zeros(nUnkEmu*(sizeEmu));
paramsODE.A1 = zeros(nUnkEmu*(sizeEmu));
for k = 2 : sizeEmu+1
    locRow = (k-2)*nUnkEmu+(1:nUnkEmu);
    locCol = (k-2)*nUnkEmu+(1:nUnkEmu);
    paramsODE.A0(locRow,locCol) = paramsODE.A0(locRow,locCol)+params.A0;
    paramsODE.A1(locRow,locCol) = paramsODE.A1(locRow,locCol)+params.A1;
end
paramsODE.B0 = zeros(nUnkEmu*(sizeEmu),nKnownEmu*(sizeEmu));
paramsODE.B1 = zeros(nUnkEmu*(sizeEmu),nKnownEmu*(sizeEmu));
for k = 2 : sizeEmu+1
    locRow = (k-2)*nUnkEmu+(1:nUnkEmu);
    locCol = (k-2)*nKnownEmu+(1:nKnownEmu);
    paramsODE.B0(locRow,locCol) = paramsODE.B0(locRow,locCol)+params.B0;
    paramsODE.B1(locRow,locCol) = paramsODE.B1(locRow,locCol)+params.B1;
end
paramsODE.Sv0 = diag(repmat(params.Sv0,sizeEmu, 1));
paramsODE.Sv1 = diag(repmat(params.Sv1,sizeEmu, 1));
paramsODE.tmpPCoefConcs = repmat(params.tmpPCoefConcs, sizeEmu,1);
tmpParamsY =reshape(params.Y(:,2:end), numel(params.Y(:,2:end)),1);
paramsODE.Y = vertcat(tmpParamsY{:});
paramsODE.t = params.t;


options=optionsODE.final;

[sol.x, sol.y]= odeFunc(@(t,X) myode(t,X, paramsODE, sizeEmu), querryTime, initX, options);

output.emu = sol.y';
output.time = sol.x';

end

function dVars = myode(t, X, paramsODE, sizeEmu)
timeDiff = (paramsODE.t-paramsODE.t(1))-(t-paramsODE.t(1))*(1-10^-6);
idt1 = find(timeDiff>0, 1)-1;
idt2 = idt1+1;
t1 =  paramsODE.t(idt1);
t2 =  paramsODE.t(idt2);

C = paramsODE.tmpPCoefConcs * [1;t;t^2];
Y1 = paramsODE.Y(:,idt1);
Y2 = paramsODE.Y(:,idt2);
Y = (Y2-Y1)./(t2-t1).* (t-t1) + Y1;

dVars0 = (paramsODE.A0*X + paramsODE.B0 * Y - paramsODE.Sv0* X) ./ C;
dVars1 = (paramsODE.A1*X + paramsODE.B1 * Y - paramsODE.Sv1* X)*t ./ C;
dVars = dVars0 + dVars1;
end


%% Dynamics of metabolites outside the network
function poolEmuDynamics = makePoolEmuDynamics(model, querryTime, emu, sizeEmu, expData, isCorrectNatAbundanceSim)
emuNetwork = model.emuNetwork(sizeEmu);
nat13CRatio = expData.nat13CRatio;

idEmu = strcmp(emu ,emuNetwork.arrayMakeB.emu);
idMet =emuNetwork.arrayMakeB.idMet(idEmu);
poolEmuDynamics = cell(1,sizeEmu+1);
if isempty(idMet)
    error('error in solveEMUDynamics.m')
end
switch model.metInfo.massType(idMet)
    case 1 % U-13C
        for k = 1 : sizeEmu+1
            poolEmuDynamics{k} = zeros(1,length(querryTime)) + ...
                nchoosekScalar(sizeEmu,k-1)*0.99^(k-1)*0.01^(sizeEmu-k+1); % 99% 13C atom purity
        end
        
    case 2 % U-12C
        if isCorrectNatAbundanceSim
            for k = 1 : sizeEmu+1
                poolEmuDynamics{k} = zeros(1,length(querryTime)) + ...
                    nchoosekScalar(sizeEmu,k-1)* nat13CRatio(2)^(k-1)* nat13CRatio(1)^(sizeEmu-k+1);
            end
        else
            for k = 1 : sizeEmu+1
                poolEmuDynamics{k} = zeros(1,length(querryTime)) + ...
                    nchoosekScalar(sizeEmu,k-1)* 0^(k-1)* 1^(sizeEmu-k+1);
            end            
        end
        
        
    case 4  %U-13C or U-12C, the ratio is accoding to input
        if expData.time(1)==0
            tmpExpTime = [expData.time];
            tmpMDV = [expData.MDVs{idMet}];
        else
            tmpExpTime = [0,expData.time];
            initMDV = zeros(sizeEmu+1,1);
            initMDV(1) = 1;
            tmpMDV = [initMDV,expData.MDVs{idMet}];
        end
        tmpMDV(isnan(tmpMDV)) = 0;
        tmpMDV(2:end-1,:) = 0;
        tmpMDV = tmpMDV ./ repmat(sum(tmpMDV,1),size(tmpMDV,1),1);        
        for k = 1 : sizeEmu+1
            switch k
                case 1
                    poolEmuDynamics{k}= interp1(tmpExpTime, tmpMDV(1,:), querryTime);
                case sizeEmu+1
                    poolEmuDynamics{k}= interp1(tmpExpTime, tmpMDV(end,:), querryTime);
                otherwise
                    poolEmuDynamics{k}= interp1(tmpExpTime, zeros(1,length(tmpExpTime)), querryTime);
            end
        end
        
end

end

function b=nchoosekScalar(n,k)

b = prod(1:n)/prod(1:k)/prod(1:n-k);

end

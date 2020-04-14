%% Draw time course of estimated fluxes
function figH=drawFigFluxTimeCourse(model, expData, solMFA, optionsMFA, options)

options.nPlotSol = min(1);
options.isRecalcErrorbar = false;
options.typePlotErrorbar = 2;
options.typeDrawCIArea = 2;
options.isPlotAllFitting = true;
options.typeFluxAll = 2;


if ~isfield(options, 'legendList')
    options.legendList = [];
end
options.plotColor = {'r'};

options.fontSize.axis = 8;
options.fontSize.title = 12;
options.fontSize.legend=8;
options.capSize=3;
options.markerSize=3;


nFig=0;

nFig=nFig+1;
figH(nFig) = figure('Color', 'w');

for i = 1 : length(options.typeFluxAll)
    typeFlux = options.typeFluxAll(i); % flux
    tmpDrawFluxTimeCourse(...
        model, expData, solMFA(1), [], optionsMFA, options, figH, typeFlux)
end

end


function tmpDrawFluxTimeCourse(...
    model, expData, solMFA, solMFARsmpl, optionsMFA, options, figH, typeFlux)
nNetRxns = optionsMFA(1).varSet.nNetRxns;
nRxns = optionsMFA(1).varSet.nRxns;
idNetRxns = optionsMFA(1).varSet.idNetRxns;

idTypeFlux = find(typeFlux==options.typeFluxAll);
fontSize=options.fontSize;


switch typeFlux
    case 1 % fllux
        nPlotRxns = nRxns;
        fieldNameFlux = 'switchTimeFluxes';
        fieldNameFluxSD = 'switchTimeFluxSDs';
        fieldNameFluxSimCIUb = 'fluxSimCIUb';
        fieldNameFluxSimCILb = 'fluxSimCILb';
        fieldNameFluxesCI= 'fluxesCI';
        fieldNameRxnNames = 'rxnNames';
        idRxnList = 1:nRxns;
        if isfield(solMFA, 'isRxnDepNSwitchTimes') && solMFA.isRxnDepNSwitchTimes
            fieldNameMatNSwitchTimesRxns = 'matNSwitchTimesRxns';
        end
    case 2 % netFlux
        nPlotRxns = nNetRxns;
        fieldNameFlux = 'netFluxes';
        fieldNameFluxSD = 'netFluxSDs';
        fieldNameFluxSimCIUb = 'netFluxSimCIUb';
        fieldNameFluxSimCILb = 'netFluxSimCILb';
        fieldNameFluxesCI= 'netFluxesCI';
        fieldNameRxnNames = 'netRxnNames';
        idRxnList = idNetRxns;
        if isfield(solMFA, 'isRxnDepNSwitchTimes') && solMFA.isRxnDepNSwitchTimes
            fieldNameMatNSwitchTimesRxns = 'matNSwitchTimesNetRxns';
        end
end

nSpRowFlux = ceil(nPlotRxns^(4/9));
nSpColFlux = ceil(nPlotRxns^(5/9));

plotCol = jet(round(length(solMFARsmpl)*1.1));
plotCol = plotCol(length(solMFARsmpl):-1:1, :);

if options.isPlotAllFitting
    nPlotSol = options.nPlotSol; 
else
    nPlotSol = 0;
end
idPlotSol = 1:nPlotSol;

figure(figH)
for r = 1 : nPlotRxns
    switch length(options.typeFluxAll)
        case 1
            nSpRow = nSpRowFlux;
            nSpCol = nSpColFlux;
            idSp = r;
        case 2
            nSpRow = nSpColFlux;
            nSpCol = nSpRowFlux*2;
            switch typeFlux
                case 1
                    idSp = (ceil(r/nSpRowFlux)-1)*nSpRowFlux+r;
                case 2
                    idSp = (ceil(r/nSpRowFlux))*nSpRowFlux+r;
            end
    end
    subplot(nSpRow,nSpCol,idSp)
    set(gca, 'FontSize',fontSize.axis)
    hold on
    
    for i = 1
        if i == 1
            tmpSol = solMFA;
        else
            tmpSol = solMFARsmpl(idPlotSol(i-1));
        end
        if i <= length(options.plotColor{idTypeFlux})
            tmpPlotColor = options.plotColor{idTypeFlux}(i);
            lineWidth = 2;
        else
            tmpPlotColor = 'k';
            lineWidth = 0.5;
        end
        
        if isfield(tmpSol, 'switchTimes')
            tmpFullSwitchTimes = [0, tmpSol.switchTimes', expData(1).time(end)];
        else
            tmpFullSwitchTimes = [];
        end
        
        %% Draw fluxes
        if isempty(tmpFullSwitchTimes) 
            switch typeFlux
                case 1
                    plot(tmpSol.timeFluxSim, tmpSol.fluxSim(r,:), ...
                        tmpPlotColor, 'LineWidth', lineWidth, 'LineStyle', '-')
                case 2
                    plot(tmpSol.timeFluxSim, tmpSol.netFluxSim(r,:), ...
                        tmpPlotColor, 'LineWidth', lineWidth, 'LineStyle', '-')
            end
        else
            plot(tmpFullSwitchTimes, tmpSol.(fieldNameFlux)(r,:), ...
                tmpPlotColor, 'LineWidth', lineWidth, 'LineStyle', '-')
        end

    end
                
    ylimData = ylim;
    if ylimData(1) >= 0
        ylim([0, ylimData(2)])
    end
    if typeFlux==1 || model.rxnInfo.revSets(idRxnList(r))==0 
        ylim([0, ylimData(2)])
    end
    ylimData = ylim;
    tmpUbFlux = model.rxnInfo.ubFluxes(idRxnList(r));
    if ylimData(2) >= tmpUbFlux
        ylim([ylimData(1),tmpUbFlux])
    end
    ylimData = ylim;
    if ylimData(1) <= -tmpUbFlux
        ylim([-tmpUbFlux,tmpUbFlux])
    end

    titleStr =model.(fieldNameRxnNames){r};
    titleStr = strrep(titleStr, '_f', '');
    titleStr = strrep(titleStr, '_', '\_');
    title(titleStr,'FontSize',fontSize.title);
    
    xlim([-5,expData(1).time(end)+5])
    set(gca, 'XTick', [0:20:expData(1).time(end)])
end

end



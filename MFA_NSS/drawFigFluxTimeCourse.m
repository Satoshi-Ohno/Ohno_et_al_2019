function figH=drawFigFluxTimeCourse(model, expData, solMFA, solMFARsmpl, optionsMFA, options)

% field2var(optionsMFA.varSet)
% field2var(solMFA)

% model = recoverRxnMetInfo(model);

if ~isfield(options, 'legendList')
    options.legendList = [];
end
if ~isfield(options, 'plotColor')
    options.plotColor = cell(1, length(options.typeFluxAll));
    for i = 1 : length(options.typeFluxAll)
        switch i
            case 1
                options.plotColor{i} = 'b';
            case 2
                options.plotColor{i} = 'r';
            otherwise
                options.plotColor{i} = 'g';
        end
    end
end
% field2var(options)

% options.fontSize.legend=options.fontSize.title;

options.fontSize.axis = 8;
options.fontSize.title = 12;
options.fontSize.legend=8;
options.capSize=3;
options.markerSize=3;


nFig=0;

%% FigureÇÃèÄîı
nFig=nFig+1;
figH(nFig) = figure('Color', 'w');
% figH(nFig) = myfigA4('V');
% Vratio = 0.5;
% Hratio = 1;
% pos = get(figH(nFig), 'Position');
% pos(2) = pos(2)+pos(4)*(1-Vratio);
% pos(4) = pos(4)*(Vratio);
% pos(3) = pos(3)*(Hratio);
% set(figH(nFig), 'Position', pos);



% figH = myfigA4;
% figPos = get(gcf, 'Position');
% figPos(1) = figPos(1)+figPos(3)/2;
% set(gcf, 'PaperPositionMode', 'auto')
% set(gcf, 'Position',figPos)


%% ÉtÉâÉbÉNÉXÇÃê}é¶
for i = 1 : length(options.typeFluxAll)
    typeFlux = options.typeFluxAll(i); % flux
    tmpDrawFluxTimeCourse170901(...
        model, expData, solMFA(1), [], optionsMFA, options, figH, typeFlux)
end

end


function tmpDrawFluxTimeCourse170901(...
    model, expData, solMFA, solMFARsmpl, optionsMFA, options, figH, typeFlux)
%% èÄîı
% field2var(optionsMFA.varSet)
% field2var(sol)
nNetRxns = optionsMFA(1).varSet.nNetRxns;
nRxns = optionsMFA(1).varSet.nRxns;
idNetRxns = optionsMFA(1).varSet.idNetRxns;
% idRxns = optionsMFA(1).varSet.idRxns;

% nKnots = optionsMFA(1).nKnots
% field2var(options)

idTypeFlux = find(typeFlux==options.typeFluxAll);
fontSize=options.fontSize;

% axisFontSize = 8;
% titleFontSize = 10;
% labelFontSize = 10;
% legendFontSize = 8;
% markerSize = 5;

%% ê}é¶èåèÇ…ÇÊÇÈê›íË
switch typeFlux
    case 1 % fllux
        nPlotRxns = nRxns;
        fieldNameFlux = 'knotFluxes';
        fieldNameFluxSD = 'knotFluxSDs';
        fieldNameFluxSimCIUb = 'fluxSimCIUb';
        fieldNameFluxSimCILb = 'fluxSimCILb';
        fieldNameFluxesCI= 'fluxesCI';
        fieldNameRxnNames = 'rxnNames';
        idRxnList = 1:nRxns;
        if isfield(solMFA, 'isRxnSpecificNKnots') && solMFA.isRxnSpecificNKnots
            fieldNameMatNKnotsRxns = 'matNKnotsRxns';
%             matNKnotsRxns = optionsMFA.matNKnotsRxns;
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
        if isfield(solMFA, 'isRxnSpecificNKnots') && solMFA.isRxnSpecificNKnots
            fieldNameMatNKnotsRxns = 'matNKnotsNetRxns';
%             matNKnotsRxns = optionsMFA.matNKnotsNetRxns;
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

%% ê}é¶
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
        otherwise
            error('typeFluxAllÇ™3à»è„ÇÃèÍçáÇÕñ¢ëŒâû') 
    end
    subplot(nSpRow,nSpCol,idSp)
    set(gca, 'FontSize',fontSize.axis)
%     set(gca, 'FontSize',fontSize.axis*0.8)
    hold on
    
    for i = 1
%     for i = nPlotSol+1:-1:1
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
        
        if isfield(tmpSol, 'knots')
            tmpFullKnots = [0, tmpSol.knots', expData(1).time(end)];
        else
            tmpFullKnots = [];
        end
        
        %% êMóäãÊä‘ÇÃê}é¶
        isDrawCIArea = false;
        switch options.typeDrawCIArea
            case 0 %ê}é¶ÇµÇ»Ç¢
            case 1 %ç≈ìKâÇÃÇ›ê}é¶
                if i== 1
                    isDrawCIArea = true;
                end
            case 2 % ëSÇƒÇ…Ç¬Ç¢Çƒê}é¶
                isDrawCIArea = true;
        end
        if ~isfield(solMFA, 'knotFluxSDs')
           isDrawCIArea = false;
        end
        if isDrawCIArea
            fillX = [tmpSol.timeFluxSim, tmpSol.timeFluxSim(end:-1:1)];
            fillY = [...
                tmpSol.(fieldNameFluxSimCIUb)(r,:),...
                tmpSol.(fieldNameFluxSimCILb)(r,end:-1:1)];
            
            fillH = fill(fillX, fillY, tmpPlotColor);
            set(fillH, 'FaceAlpha', 0.3, 'EdgeColor', 'none')
        end

        %% ÉGÉâÅ[ÉoÅ[ÇÃê}é¶
        isPlotErrorbar = false;
        switch options.typePlotErrorbar
            case 0 %ê}é¶ÇµÇ»Ç¢
            case 1 %ç≈ìKâÇÃÇ›ê}é¶
                if i== 1
                    isPlotErrorbar = true;
                end
            case 2 % ëSÇƒÇ…Ç¬Ç¢Çƒê}é¶
                isPlotErrorbar = true;
        end
        if ~isfield(solMFA, 'knotFluxSDs')
           isPlotErrorbar = false;
        end
        if isempty(tmpFullKnots)
            isPlotErrorbar = false;
        end
        if isPlotErrorbar
            loc = ismember(tmpSol.timeFluxSim, tmpFullKnots);
            halfCI = (tmpSol.(fieldNameFluxSimCIUb)(r,loc)+tmpSol.(fieldNameFluxSimCILb)(r,loc))/2;
            lenCI = (tmpSol.(fieldNameFluxSimCIUb)(r,loc)-tmpSol.(fieldNameFluxSimCILb)(r,loc))/2;
            if isfield(tmpSol, 'isRxnSpecificNKnots') && tmpSol.isRxnSpecificNKnots
                isPlotErrorbarKnots = optionsMFA(i).(fieldNameMatNKnotsRxns)(r,:);
%                 isPlotErrorbarKnots = matNKnotsRxns(r,:);
            else
                isPlotErrorbarKnots = true(1,length(tmpFullKnots));
            end
            errorbar(tmpFullKnots(isPlotErrorbarKnots), halfCI(isPlotErrorbarKnots) , lenCI(isPlotErrorbarKnots), ...
                tmpPlotColor, 'LineWidth', lineWidth, 'LineStyle', 'none')
        end
        
        %% ÉtÉâÉbÉNÉXÇÃê}é¶
        if isempty(tmpFullKnots) % LakeEeóp
            switch typeFlux
                case 1
                    plot(tmpSol.timeFluxSim, tmpSol.fluxSim(r,:), ...
                        tmpPlotColor, 'LineWidth', lineWidth, 'LineStyle', '-')
                case 2
                    plot(tmpSol.timeFluxSim, tmpSol.netFluxSim(r,:), ...
                        tmpPlotColor, 'LineWidth', lineWidth, 'LineStyle', '-')
            end
        else
            plot(tmpFullKnots, tmpSol.(fieldNameFlux)(r,:), ...
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
    tmpUbFlux = optionsMFA(1).ub.knotFluxes;
%     tmpUbFlux = max([optionsMFA(1).ub.knotFluxes, optionsMFA(1).constr.fluxes.ub]);
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
%     title(titleStr,'FontSize',titleFontSize*0.8);
    
    xlim([-5,expData(1).time(end)+5])
    set(gca, 'XTick', [0:20:expData(1).time(end)])
end

%% legend
% if isempty(options.legendList)
%     return
% end
% 
% subplot(nSpRow,nSpCol,idSp+1)
% hold on
% % set(gca, 'FontSize',legendFontSize*1.2, 'Visible', 'off')
% for i = 1:nPlotSol+1
%     if i <= length(options.plotColor{idTypeFlux})
%         tmpPlotColor = options.plotColor{idTypeFlux}(i);
%         lineWidth = 2;
%     else
%         tmpPlotColor = 'k';
%         lineWidth = 0.5;
%     end
%     plot(1:2, 1:2, tmpPlotColor, 'LineWidth', lineWidth, 'LineStyle', '-')
% %     plot(1:10, 1:10, tmpPlotColor, 'LineWidth', lineWidth, 'LineStyle', '-', 'Visible', 'off')
% end
% xlim([-1e1,1e1])
% ylim([-1e1,1e1])
% set(gca, 'Visible', 'off', 'FontSize', fontSize.axis)
% legend(options.legendList,'Location','SouthWest', 'FontSize',fontSize.legend)

end



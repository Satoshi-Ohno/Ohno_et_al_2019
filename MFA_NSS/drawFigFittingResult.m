%% Draw fitting results for metabolite concentrations and mass isotopomer fractions
function     figH = drawFigFittingResult(model, expData, solMFA, optionsMFA, options)

nExpType = length(solMFA);
nPlotMet = length(options.plotMet);
options.plotMet = model.mets(optionsMFA(1).varSet.idNonPoolMets)';

[~, idPlotMet]= ismember(options.plotMet, model.mets);
if any(idPlotMet==0)
    error("Any of options.plotMet are not inculded in the model.mets")
end

options.fontSize.axis = 8;
options.fontSize.title = 12;
options.fontSize.legend=8;
options.capSize=3;
options.markerSize=3;

options.colConcs = [0,0,0];
options.colMDV = jet(7);

idCol = 0;
nFig = 0;
for e = 1 : nExpType
    for i = 1 : 2 
        nFig = nFig + 1;
        figH(nFig) = figure('Color', 'w');
        nSpRow = ceil(nPlotMet^(4/9));
        nSpCol = ceil(nPlotMet^(5/9));
        
        for mm = 1 : nPlotMet+1
            idCol = mm;
            if mm <= nPlotMet
                m = idPlotMet(mm);
                options.isLegend=false;
            else
                m = idPlotMet(end);
                options.isLegend=true;
            end
            
            options.idExpType = e;
            options.idMet = m;
            if i==1 && options.isLegend
                continue
            end
            switch i                
                case 1 % metabolit concentrations
                    options.plotType = 'conc';
                    
                case 2 % mass isitopomer fractions
                    options.plotType  = 'MDV';
            end
            idSp = mm;
            spH(e,mm,i) = subplot(nSpRow, nSpCol, idSp);
            hold on
            drawH = drawExpSim(model,expData, solMFA, optionsMFA, options);
            set(gca, 'Color', 'none')
        end
    end
end

for mm = 1: nPlotMet
    maxYLim = 0;
    for e = 1 : nExpType
        tmpYLim = get(spH(e,mm,1), 'YLim');
        maxYLim = max([maxYLim, tmpYLim(2)]);
    end
    for e = 1 : nExpType
        set(spH(e,mm,1), 'YLim', [0,maxYLim]);
    end    
end



end

%% Draw concentrations and mass isotopomer fractions for each metabolite
function drawH = drawExpSim(model,expData, solMFA, optionsMFA, options)

nDraw = 0;
e = options.idExpType;
m = options.idMet;
switch options.plotType
    %% metabolite concentrations
    case {'conc'}
        xExp = expData(e).time;
        yExp = expData(e).concs(m,:);
        yExpErrbar = expData(e).concsSD(m,:);
        
        xSim = solMFA(e).timeSim;
                
        ySim = solMFA(e).concs(m,:);
        nDraw = nDraw +1;
        drawH(nDraw)= plot(xSim, ySim,  '-', 'LineWidth', 1, 'Color', options.colConcs);
        
        if optionsMFA(1).varSet.isEvalConcMets(m)
            nDraw = nDraw +1;
            drawH(nDraw) = errorbar(xExp, yExp, yExpErrbar, 'o', 'Color', options.colConcs);
            set(drawH(nDraw), ...
                'CapSize', options.capSize, ...
                'LineStyle', 'none',...
                'MarkerSize', options.markerSize)
        end
        
        ySimExpTime = selectCorresData(ySim, xSim, xExp,2);
        loc = ~isnan(yExp);
        tmpCorrCoef = corrcoef([yExp(loc)', ySimExpTime(loc)']);
        Rsq = tmpCorrCoef(1,2).^2;
        
        ylimVal = ylim;
        ylim([0,ylimVal(2)]);
        
    %% mass isotopomer fractions
    case {'MDV'}
        xExp = expData(e).time;
        yExp = expData(e).MDVs{m};
        yExpErrbar = expData(e).MDVsSD{m};
        if options.isLegend
            for j = 1 : size(options.colMDV,1)
                nDraw = nDraw+1;
                drawH(nDraw) = errorbar(xExp, yExp(1,:), yExpErrbar(1,:), ['o']);
                set(drawH(nDraw), ...
                    'Color', options.colMDV(j,:),...
                    'CapSize', options.capSize, ...
                    'LineStyle', '-',...
                    'MarkerSize', options.markerSize)
            end
            
        else
            xSim = solMFA(e).timeSim;
            
        
            ySim = solMFA(e).MDVsSim{m};
            for j = 1 : size(ySim,1)
                nDraw = nDraw +1;
                drawH(nDraw) = plot(xSim, ySim(j,:), '-', 'LineWidth', 1, 'Color', options.colMDV(j,:));
            end
            if optionsMFA(1).varSet.isEvalMDVMets(m)
                for j = 1 : size(yExp,1)
                    nDraw = nDraw +1;
                    drawH(nDraw) = errorbar(xExp, yExp(j,:), yExpErrbar(j,:), ['o']);
                    set(drawH(nDraw), ...
                        'Color', options.colMDV(j,:),...
                        'CapSize', options.capSize, ...
                        'LineStyle', 'none',...
                        'MarkerSize', options.markerSize)
                end
            end
            
            
           ySimExpTime = selectCorresData(ySim, xSim, xExp,2);
            tmpYExp = yExp(:);
            tmpYSimExpTime = ySimExpTime(:);
            loc = ~isnan(tmpYExp);
            tmpCorrCoef = corrcoef([tmpYExp(loc), tmpYSimExpTime(loc)]);
            Rsq = tmpCorrCoef(1,2).^2;
            
            ylim([0,1.1])
            
        end
        
end


%% legend and setting
set(gca, 'FontSize', options.fontSize.axis)
if options.isLegend
    set(gca, 'Visible', 'off')
    xlim([-1e3,1e3])
    ylim([-1e3,1e3])
    legendList = cell(size(options.colMDV,1),1);
    for i = 1:size(options.colMDV,1)
        legendList{i} = ['M+' num2str(i-1)];
    end
    legendH = legend(legendList, 'Location', 'SouthWest', 'FontSize', options.fontSize.legend);
else
    xlim([0-xExp(end)*0.05,xExp(end)+xExp(end)*0.05])
    set(gca, 'XTick', 0:20:xExp(end))
    set(gca, 'TickLength', get(gca, 'TickLength')*3)
    
    titleStr = [model.mets{m} ' '];
    titleStr = strrep(titleStr, '_', '\_');
    title(titleStr, 'FontSize', options.fontSize.title)
    box off
    
end
end


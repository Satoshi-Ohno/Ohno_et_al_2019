function     figH = drawFigFittingResult(model, expData, solMFA, solMFARsmpl, optionsMFA, options)
% 代謝物濃度のフィッティング結果
% 同位体割合のフィッティング結果

nExpType = length(solMFA);
nPlotMet = length(options.plotMet);
typeFigStr = options.typeFigStr;

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


%% 図示
% figH = myfigA4('V');
% set(figH, 'Color', 'none')
if ismember({'2'}, typeFigStr)
    Vratio = 0.1;
    Hratio = 0.6;
elseif ismember({'S1'}, typeFigStr)
    Vratio = 0.4;
    Hratio = 0.6;
elseif ismember({'S2'}, typeFigStr)
    Vratio = 0.4;
    Hratio = 0.6;
end
nSpRow = 6;

idCol = 0;
% for e = options.idExpType
nFig = 0;
for e = 1 : nExpType
    for i = 1 : 2 % plotType
        nFig = nFig + 1;
%         figH(nFig) = myfigA4('V');
%         % set(figH, 'Color', 'none')
%         pos = get(figH(nFig), 'Position');
%         pos(2) = pos(2)+pos(4)*(1-Vratio);
%         pos(4) = pos(4)*(Vratio);
%         pos(3) = pos(3)*(Hratio);
%         set(figH(nFig), 'Position', pos);
        figH(nFig) = figure('Color', 'w');
        nSpRow = ceil(nPlotMet^(4/9));
        nSpCol = ceil(nPlotMet^(5/9));

%         if ismember({'3'}, typeFigStr)
%             nSpCol = nPlotMet+1;
%         else
%             nSpCol = 5;
%         end
        
        for mm = 1 : nPlotMet+1
%             if e==1 && mm == nPlotMet+1
%                 continue
%             end
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
                %% 代謝物濃度
                case 1
                    options.plotType = 'conc';
                    %% 同位体割合
                case 2
                    options.plotType  = 'MDV';
            end
            if ismember({'3'}, typeFigStr)
                idRow = [2,3,4];
                if options.isLegend
                    idSp = sub2ind([nSpCol, nSpRow], [idCol,idCol,idCol],[idRow]);
                else
                    idSp = sub2ind([nSpCol, nSpRow], [idCol,idCol,idCol],idRow);
                end
            else
                idSp = mm;
            end
%             if e==2 && mm == nPlotMet+1
%                 keyboard;
%             end
            spH(e,mm,i) = subplot(nSpRow, nSpCol, idSp);
            hold on
            drawH = drawExpSim(model,expData, solMFA, solMFARsmpl, optionsMFA, options);
            set(gca, 'Color', 'none')
        end
    end
end

%% 軸の範囲をそろえる

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

%% 代謝物濃度・同位体割合の図示
function drawH = drawExpSim(model,expData, solMFA, solMFARsmpl, optionsMFA, options)

nDraw = 0;
e = options.idExpType;
m = options.idMet;
switch options.plotType
    case {'conc'}
        xExp = expData(e).time;
        yExp = expData(e).concs(m,:);
        yExpErrbar = expData(e).concsSD(m,:);
        
        xSim = solMFA(e).timeSim;
        
        % 信頼区間の図示
        if ~isempty(solMFARsmpl)
        [ySimCIUb, ySimCILb] =  calcSimCI(model,expData, solMFA, solMFARsmpl, options);
        nDraw = nDraw +1;
        fillX = [xSim, xSim(end:-1:1)];
        fillY = [ySimCIUb, ySimCILb(end:-1:1)];
        fillH = fill(fillX, fillY, 'k');
        tmpCol = options.colConcs;
%         tmpCol=tmpCol+0.7; % 色を白っぽくする
%         tmpCol(tmpCol>1) = 1;
        set(fillH, 'FaceColor', tmpCol)
%         set(fillH, 'EdgeColor', 'none')
        set(fillH, 'FaceAlpha', 0.3, 'EdgeColor', 'none')
        end
        
        % 推定値の図示
        ySim = solMFA(e).concs(m,:);
        nDraw = nDraw +1;
        drawH(nDraw)= plot(xSim, ySim,  '-', 'LineWidth', 1, 'Color', options.colConcs);
        
        % 実験データの図示
        if optionsMFA(1).varSet.isEvalConcMets(m)
            nDraw = nDraw +1;
            drawH(nDraw) = errorbar(xExp, yExp, yExpErrbar, 'o', 'Color', options.colConcs);
            set(drawH(nDraw), ...
                'CapSize', options.capSize, ...
                'LineStyle', 'none',...
                'MarkerSize', options.markerSize)
        end
        
        % R^2値の算出
        ySimExpTime = selectCorresData(ySim, xSim, xExp,2);
        loc = ~isnan(yExp);
        tmpCorrCoef = corrcoef([yExp(loc)', ySimExpTime(loc)']);
        Rsq = tmpCorrCoef(1,2).^2;
        
        ylimVal = ylim;
        ylim([0,ylimVal(2)]);
        
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
            
            % 信頼区間の図示
        if ~isempty(solMFARsmpl)
            [ySimCIUb, ySimCILb] =  calcSimCI(model,expData, solMFA, solMFARsmpl, options);
            nDraw = nDraw +1;
            fillX = [xSim, xSim(end:-1:1)];
            for j = 1 : size(ySimCIUb,1)
                fillY = [ySimCIUb(j,:), ySimCILb(j,end:-1:1)];
                fillH = fill(fillX, fillY, 'k');
                tmpCol = options.colMDV(j,:);
%                 tmpCol=tmpCol+0.7; % 色を白っぽくする
%                 tmpCol(tmpCol>1) = 1;
                set(fillH, 'FaceColor', tmpCol)
%                 set(fillH, 'EdgeColor', 'none')
                set(fillH, 'FaceAlpha', 0.3, 'EdgeColor', 'none')
            end
        end
        
            ySim = solMFA(e).MDVsSim{m};
            for j = 1 : size(ySim,1)
                nDraw = nDraw +1;
                drawH(nDraw) = plot(xSim, ySim(j,:), '-', 'LineWidth', 1, 'Color', options.colMDV(j,:));
            end
            % 実験値の図示
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
            
            % R^2値の算出
            ySimExpTime = selectCorresData(ySim, xSim, xExp,2);
            tmpYExp = yExp(:);
            tmpYSimExpTime = ySimExpTime(:);
            loc = ~isnan(tmpYExp);
            tmpCorrCoef = corrcoef([tmpYExp(loc), tmpYSimExpTime(loc)]);
            Rsq = tmpCorrCoef(1,2).^2;
            
            ylim([0,1.1])
            
        end
        
end

%% 余白を削る
% https://jp.mathworks.com/help/matlab/creating_plots/save-figure-with-minimal-white-space.html#butpbwd-1
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];

%% legendと設定
% set(gca, 'FontName', options.fontName)
set(gca, 'FontSize', options.fontSize.axis)
if options.isLegend
    set(gca, 'Visible', 'off')
    xlim([-1e3,1e3])
    ylim([-1e3,1e3])
%     xlim([-1e6,1e6])
%     ylim([-1e6,1e6])
    legendList = cell(size(options.colMDV,1),1);
    for i = 1:size(options.colMDV,1)
        legendList{i} = ['M+' num2str(i-1)];
    end
    legendH = legend(legendList, 'Location', 'SouthWest', 'FontSize', options.fontSize.legend);
%     legendH = legend(legendList, 'Location', 'SouthWest', 'FontSize', options.fontSize.legend, 'FontName',options.fontName);
    if ismember({'3'}, options.typeFigStr)
        pos = get(legendH, 'Position');
        %     pos(1)=0;
        pos(2)=0.1;
        set(legendH, 'Position', pos)
    end
else
    xlim([0-xExp(end)*0.05,xExp(end)+xExp(end)*0.05])
    set(gca, 'XTick', 0:20:xExp(end))
    set(gca, 'TickLength', get(gca, 'TickLength')*3)
    
    titleStr = [model.mets{m} ' '];
    titleStr = strrep(titleStr, '_', '\_');
%     titleStr = [titleStr '(' num2str(Rsq,'%1.2f') ')'];
%     titleStr = char({titleStr, ['\fontsize{' num2str(options.fontSize.title*0.8) '} R^2=', num2str(Rsq,'%1.2f')]});
%     titleStr = strjust(titleStr);
%     title(titleStr)
    title(titleStr, 'FontSize', options.fontSize.title)
%     title(titleStr, 'FontSize', options.fontSize.title, 'FontName', options.fontName)
    box off
    
%     if isnan(Rsq) || Rsq<=0.5
%         keyboard;
%     end
end
end

%% 信頼区間の算出
function   [ySimCIUb, ySimCILb] =  calcSimCI(model,expData, solMFA, solMFARsmpl, options)
nDraw = 0;
e = options.idExpType;
m = options.idMet;

alphaCI = options.alphaCI;
N = size(solMFARsmpl,2);
nSmplCI = N*alphaCI;
idSolCalcCI = (1+(N-nSmplCI)/2) : ((N-nSmplCI)/2 + nSmplCI);
idSolCalcCI = ceil(idSolCalcCI);

switch options.plotType
    case {'conc'}
        ySimRsmpl = cat(3, solMFARsmpl(e,:).concs);
        ySimRsmpl = ySimRsmpl(m,:,:);
        ySimRsmpl = sort(ySimRsmpl,3);
        ySimCIUb = ySimRsmpl(:,:,idSolCalcCI(end));
        ySimCILb = ySimRsmpl(:,:,idSolCalcCI(1));
        
    case {'MDV'}
        [nRow, nCol] = size(solMFARsmpl(e,1).MDVsSim{m});
        ySimRsmpl = zeros(nRow,nCol,N);
        for r = 1 : N
            ySimRsmpl(:,:,r) = solMFARsmpl(e,r).MDVsSim{m};
        end
        ySimRsmpl = sort(ySimRsmpl,3);
        ySimCIUb = ySimRsmpl(:,:,idSolCalcCI(end));
        ySimCILb = ySimRsmpl(:,:,idSolCalcCI(1));
end

end

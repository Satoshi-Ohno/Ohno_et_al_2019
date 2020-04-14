%% Draw scatter plots of input and estimated fluxes
function figH = drawScatterInputVsEstimatedFlux(solKAMinAIC, solKANoReg, optionsKM)

timeKinetic = optionsKM.timeKinetic;
nKineticTime = length(timeKinetic);

fontSize.axis = 10;
fontSize.title = 15;
fontSize.legend=10;
fontSize.label=15;

nFig = 0;

for i = 1 : 2
    nFig = nFig+1;
    figH(nFig) = figure('Color', 'w');
    if i==1
        titleStr = 'Model with only substrates and products';
    else
        titleStr = 'Model with minimized AIC';
    end
    
    if i == 1
        tmpSol = solKANoReg;
    else
        tmpSol = solKAMinAIC;
    end
    normExpData = optionsKM.normExpData;
    locFlux = strcmp('flux', {normExpData.regulatorType});
    fluxMFA = colVec(normExpData(locFlux).level)';
    fluxEst = colVec(tmpSol.fluxEst)';
    
    fluxEst=fluxEst(fluxMFA~=0);
    fluxMFA=fluxMFA(fluxMFA~=0);
    
    hold on
    
    if length(fluxMFA) <= nKineticTime*2
        xIns = log(abs(fluxMFA(1:nKineticTime)));
        yIns = log(abs(fluxEst(1:nKineticTime)));
        xCtrl = log(abs(fluxMFA(nKineticTime+(1:nKineticTime))));
        yCtrl = log(abs(fluxEst(nKineticTime+(1:nKineticTime))));
    else
        fluxMFA = fluxMFA(1:nKineticTime*2) - fluxMFA(nKineticTime*2+1:end);
        fluxEst = fluxEst(1:nKineticTime*2) - fluxEst(nKineticTime*2+1:end);
        xIns = log(abs(fluxMFA(1:nKineticTime)));
        yIns = log(abs(fluxEst(1:nKineticTime)));
        xCtrl = log(abs(fluxMFA(nKineticTime+(1:nKineticTime))));
        yCtrl = log(abs(fluxEst(nKineticTime+(1:nKineticTime))));
    end
    plot(xIns, yIns, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
    plot(xCtrl, yCtrl, 'b^', 'MarkerSize', 5, 'MarkerFaceColor', 'b');
    plot([-100,100],[-100,100], 'k-')
    
    xRange = [min(log(abs(fluxMFA))), max(log(abs(fluxMFA)))];
    yRange = [min(log(abs(fluxEst))),max(log(abs(fluxEst)))];
    xyRange = [min([xRange(1),yRange(1)]), max([xRange(2),yRange(2)])];
    xlim(xyRange)
    ylim(xyRange)
    
    % reguression line
    X = [ones(length(fluxMFA),1), log(abs(fluxMFA'))];
    Y = log(abs(fluxEst'));
    regCoef = X\Y;
    Ycalc = X*regCoef;
    plot(log(abs(fluxMFA)), Ycalc, '-', 'Color', [0.8,0,0.8]);
    
    
    set(gca, 'FontSize', fontSize.axis)
    xlabel(' log(input flux)', 'FontSize', fontSize.label)
    ylabel('log(kinetic flux)', 'FontSize', fontSize.label)
    
    tmp = corrcoef(X(:,2), Y);
    Rsq = tmp(1,2).^2;
    
    if isnan(Rsq)
        Rsq = 0;
    end
    textStr =  char(...
        [' R^2 = ' num2str(Rsq,'%1.2f')], ...
        ['AIC = ' num2str(tmpSol.AIC, '%1.1f')]);
    textH=text(xyRange(1) + 0.1*range(xyRange), xyRange(2)- 0.2*range(xyRange),...
        textStr);
    set(textH, 'FontSize',fontSize.legend)
    title(titleStr, 'FontSize', fontSize.title);
    
    
    
    legendList = {
        'Ins', ...
        'Ctrl',...
        'y=x', 'regression line'};
    legend(legendList, 'Location','SouthEast', 'FontSize',fontSize.legend)
end


end


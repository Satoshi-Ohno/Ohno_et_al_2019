%% Draw bar graph of time-averaged regulation coefficient
function figH = drawRegCoef(regCoefFxn, optionsKA);

fxnRegulatorNames  = regCoefFxn(1).fxnRegulatorNames;
legendList = regCoefFxn(1).legendList;

fontSize.axis = 10;
fontSize.title = 15;
fontSize.legend=10;
fontSize.label=15;

nFig = 0;

nFig = nFig+1;
figH(nFig) = figure('Color', 'w');
titleStr = 'Time courses of regulation coefficients';

lineStyle = '-';

hold on
    
colMatFxnReg = [
    0,1,0; % green: phosphorylation
    1,0,0; % redÅFallostery
    0,0,1; % blueÅFsubstrates and products
    0.4,0.4,0.4;]; % gray: unaccounted regulators

tmpRegCoefFxn = regCoefFxn.rhoMat;

timeRho = regCoefFxn.timeRho;
for f = 1 : length(fxnRegulatorNames)
    plot(timeRho, tmpRegCoefFxn(f,:), lineStyle, 'Color',colMatFxnReg(f,:));
end
plot([0,60], [1,1], ':', 'Color', [0.6,0.6,0.6])
        
set(gca, 'FontSize', fontSize.axis)

xlabel('time (min)', 'FontSize', fontSize.label)
ylabel('\rho', 'FontSize', fontSize.label)
title(titleStr, 'FontSize', fontSize.title);
    
legend(legendList, 'Location','NorthEast', 'FontSize', fontSize.legend)

end
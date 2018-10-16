%% Draw bar graph of time-averaged regulation coefficient
function figH = drawTimeAveragedRegCoef(regCoefFxn, optionsKA)

fxnRegulatorNames  = regCoefFxn.fxnRegulatorNames;
legendList = regCoefFxn.legendList;

fontSize.axis = 10;
fontSize.title = 15;
fontSize.legend=10;
fontSize.label=15;

nFig = 0;

nFig = nFig+1;
figH(nFig) = figure('Color', 'w');

titleStr = ['Time-averaged regulation coefficient'];
colMatFxnReg = [
    0,1,0; % green: phosphorylation
    1,0,0; % redÅFallostery
    0,0,1; % blueÅFsubstrates and products
    0.4,0.4,0.4;]; % gray: unaccounted regulators

tmpRhoMatTimeAverage = regCoefFxn.rhoMatTimeAverage;

hold on
for f = 1 : length(fxnRegulatorNames)
    barH = bar(f, tmpRhoMatTimeAverage(f), 'FaceColor',colMatFxnReg(f,:), 'EdgeColor', 'k');
end
plot([0,length(fxnRegulatorNames)+1], [1,1], ':', 'Color', [0.6,0.6,0.6])


xlim([0, length(fxnRegulatorNames)+1])
set(gca, 'XTick', 1:length(fxnRegulatorNames), 'XTickLabel', legendList)
ylim([0, 1.1])
set(gca,'YTick', [0,0.5,1])

title(titleStr, 'FontSize', fontSize.title);
ylabel('time-averaged \rho', 'FontSize', fontSize.label)
set(gca, 'FontSize', fontSize.label)

end
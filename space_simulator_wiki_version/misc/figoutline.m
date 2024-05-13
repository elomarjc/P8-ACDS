FigPos = [0.7, 0.8, 0.5, 0.4];
FigWinSize = [8, 4];

% Position of figure and size [left bottom width height]
PaperPos = [0.25 7.7 FigWinSize];

set(0, 'Units', 'Pixels');
ScrSizePix = get(0,'ScreenSize');
set(0,'Units', 'Inches');
ScrSizeInc = get(0,'ScreenSize');
set(0, 'Units', 'Pixels');
Inc2Pix = ScrSizePix(3)/ScrSizeInc(3);
Xpos = (ScrSizeInc(3)-FigWinSize(1))/2;
Ypos = ScrSizeInc(4)-round(0.15*ScrSizeInc(4))-FigWinSize(2);
pos = Inc2Pix.*[Xpos, Ypos, FigWinSize(1), FigWinSize(2)];

% Plot data
figure('Position', pos);
set(gcf, 'DefaultLineLineWidth', 1.5);
set(gcf, 'DefaultAxesFontSize', 12);
set(gcf, 'DefaultAxesFontName', 'Times');
set(gcf, 'DefaultTextFontSize', 12);
set(gcf, 'DefaultTextFontName', 'Times');
set(gcf,'DefaultAxesColorOrder',[0 0 0],...
    'DefaultAxesLineStyleOrder','-|:|--|-.')
axes('Position', [FigPos(1)/FigWinSize(1), ...
    FigPos(3)/FigWinSize(2), ...
    (FigWinSize(1)-FigPos(1)-FigPos(2))/FigWinSize(1), ...
    (FigWinSize(2)-FigPos(3)-FigPos(4))/FigWinSize(2)]);

% Same size as window
%set(gcf,'PaperPositionMode','auto')

% Same size as window, but moved to top of paper
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', PaperPos);


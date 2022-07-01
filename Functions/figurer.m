function col = figurer
% Function that standardizes plots and makes them look good.

fontsize         = 12;
% boundingbox      = [45 3 501 501];
% boundingbox      = [0 0 100 100];
% boundingboxlarge = [45 3 549 549];
% viewport         = [81 81 390 390];
% viewport         = [0 0 100 100];
% paperposition    = [0 0 boundingbox(3:4)];

set(0,'DefaultFigureColor','w');
% set(0,'DefaultAxesColor','none');
set(0,'DefaultLegendBox','off');
set(0,'DefaultLegendInterpreter','latex');
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesTickLabelInterpreter','latex');
set(0,'DefaultTextFontsize',fontsize);
set(0,'DefaultAxesFontsize',fontsize);
set(0,'DefaultAxesLinewidth',0.1);
set(0,'DefaultLineLinewidth',0.2);
set(0,'DefaultHistogramLinewidth',0.1);
set(0,'DefaultStairLineWidth',0.1);
set(0,'DefaultErrorBarLineWidth',0.1);
set(groot,'defaultfigureposition',[200 200 550 550])
% set(0, 'box', 'on')
% set(0,'DefaultFigureUnits','points');
% set(0,'DefaultFigurePosition',boundingbox);
% set(0,'DefaultFigureResize','off');
% set(0,'DefaultAxesActivePositionProperty','Position');
% set(0,'DefaultAxesUnits','points');
% set(0,'DefaultAxesPosition',viewport);

% Define color vector for plotting
col     = ["-" ,"-" ,"-" ,"-" ,"-" ,"-" ,"-", ...
           "--","--","--","--","--","--","--",...
           "-*","-*","-*","-*","-*","-*","-*"];




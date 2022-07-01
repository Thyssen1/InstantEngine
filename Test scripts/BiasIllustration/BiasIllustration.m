%% Initialization
clear; close all; clc

% Add path to experimental plan MATLAB function
addpath(strcat(erase(pwd, "Test scripts\BiasIllustration"), "\Functions"))

figurer;

rng(10, 'twister')

%% 
bool_SaveFigures = false;

%% Bias-Noise graphic
N = 50;
X = zeros(N,2,4);

mu    = [0; 0.5];
sigma = [0.05; 0.20];

XTitle_Bias = ["No bias"; "Bias"];
YTitle_Noise = ["No noise"; "Noise"];

figure
for i = 1:2
    for j = 1:2
        X(:,:,2*(i-1)+j) = normrnd(mu(j), sigma(i), N, 2);
        
        subplot(2,2,2*(i-1)+j); hold all
        plot(0, 0, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r')
        scatter(X(:,1,2*(i-1)+j), X(:,2,2*(i-1)+j), 'bx')
        xlim([-1 1])
        ylim([-1 1])
        xlabel('$X_{1}$', 'interpreter', 'latex')
        ylabel({YTitle_Noise(i), '$X_{2}$'}, 'interpreter', 'latex')
        title(XTitle_Bias(j))
    end
end
FigureTitle(1) = "InstantEngine_BiasGraphic1";

%% Biased sine function
rng(10, 'twister')
x = (0:(2*pi/1000):2*pi)';

XVal{1} = pi*rand(10,1)+1/2*pi;
XVal{2} = 2*pi*rand(100,1);
XVal{3} = 2*pi*rand(1000,1);

fun = @(t) -sin(t);

for i = 1:3
        lm{i} = fitlm(XVal{i}, fun(XVal{i}));
end

YMin = min(min([sin(x) predict(lm{1},x) predict(lm{2},x) predict(lm{3},x)]));
YMax = max(max([sin(x) predict(lm{1},x) predict(lm{2},x) predict(lm{3},x)]));

figure
for i = 1:3
    subplot(3,1,i); hold all

    plot(x, fun(x), '-')
    plot(XVal{i}, fun(XVal{i}), '.')
    plot(x, predict(lm{i}, x), '--')
    xlim([min(x) max(x)])
    ylim([YMin YMax])
    xlabel('x')
    ylabel('f(x)')

    stats{i} = rs_stats(fun(x), predict(lm{i}, x));
    R2(i,1) = stats{i}.R2;
    RMSE(i,1) = stats{i}.RMSE;

    legend("True function", "Data", "Fit", 'location', 'northwest', ...
        'FontSize', 8, 'box', 'on')
end

disp([R2 RMSE])
FigureTitle(2) = "InstantEngine_BiasGraphic2";

%% 
for i = 1:length(FigureTitle)
    saveas(figure(i), strcat(pwd, "\Images\", FigureTitle(i)), 'png')
end



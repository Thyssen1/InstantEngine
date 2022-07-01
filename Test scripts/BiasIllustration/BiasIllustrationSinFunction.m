



x = 0:0.01:(2*pi);
Xs1 = lhsdesign(5,1);
Xs2 = lhsdesign(250,1);
Xs1 = Xs1*2*pi;
Xs2 = Xs2*2*pi;

% Fit linear models to data
mdl1 = fitlm(Xs1, sin(Xs1));
mdl2 = fitlm(Xs2, sin(Xs2));

figure
subplot(2,1,1); hold all
plot(x, sin(x), '-')
plot(Xs1, sin(Xs1), 'o')
plot(x, predict(mdl1, x(:)), '-')
legend("Ground truth", "5 samples",  "Fit 5 samples")

subplot(2,1,2); hold all
plot(x, sin(x), '-')
plot(Xs2, sin(Xs2), 'o')
plot(x, predict(mdl1, x(:)), '-')
plot(x, predict(mdl2, x(:)), '-')
legend("Ground truth", "25 samples", "Fit 5 samples", "Fit 250 samples")


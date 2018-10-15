%% HHGDEMO.M
% A script demonstrating the usage of HHGPermutationTest.m

% Set the sample size
n = 1000;

% Set the number of permutations
nperm = 100;

%% Linear correlation

x = randn(1000,1);
y = 10*x+1;

% Run the test
[p, t, pstat] = HHGPermutationTest(x,y, nperm);

% Visualize the results
hhgDemoPlot(p,t,pstat,nperm,'Perfect linear correlation');

%% Quadratic

x = randn(1000,1);
y = x.^2;

% Run the test
[p, t, pstat] = HHGPermutationTest(x,y, nperm); 

% Visualize the results
hhgDemoPlot(p,t,pstat,nperm,'Quadratic');

%% Exponential

x = randn(1000,1);
y = exp(x);

% Run the test
[p, t, pstat] = HHGPermutationTest(x,y, nperm); 

% Visualize the results
hhgDemoPlot(p,t,pstat,nperm,'Exponential');

%% Uncorrelated Gaussian variables

x = randn(1000,1);
y = randn(1000,1);

% Run the test
[p, t, pstat] = HHGPermutationTest(x,y, nperm); 

% Visualize the results
hhgDemoPlot(p,t,pstat,nperm,'Uncorrelated Gaussian variables');

%% Utility function to plot results

function hhgDemoPlot(p,t,pstat,nperm,titleStr)
figure; hold on;
plot(pstat, '-k','linewidth',2);
plot(xlim, [t t], '-r','linewidth',2);
legend({'test statistic under permutation','empirical value'});
xlabel('permutation');
ylabel('HHG test statistic');
if p == 0
    title(sprintf('%s: p < 1/%d', titleStr, nperm));
else
    title(sprintf('%s: p = %f', titleStr, p));
end
end
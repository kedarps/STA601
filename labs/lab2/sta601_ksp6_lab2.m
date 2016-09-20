%% STA 601: Lab 2
% Author: Kedar S Prabhudesai
% Created on: 09/10/2013

close all;clear all;
%% Part 2
% Parameters for T-Distribution
nu = 1;
% Number of samples to draw
nSamples = 10000;
% Distribution for T2
T2Dist = makedist('Gamma','a',nu/2,'b',2/nu);
% Draw 10,000 samples from T2 Dist
T2Rand = T2Dist.random(1,nSamples);
XGivenT2Rand = zeros(1,nSamples);

% For the Given T2 sample, draw corresponding X|T2 sample
for iSample = 1:nSamples
    XGivenT2Dist = makedist('Normal','mu',0,'sigma',sqrt(1/T2Rand(iSample))); 
    XGivenT2Rand(iSample) = XGivenT2Dist.random();
end
% Plot Histogram of Values
XHistFig = figure;
hist(gca,XGivenT2Rand,10000);
title('Histogram of X|\tau^2','Fontsize',14);
xlabel('Values','Fontsize',14);
ylabel('Count','Fontsize',14);
matlab2tikz(XHistFig,'XGivenT2Hist.tikz');

%% Part 3
% Marginal Distribution of X
XDist = makedist('tLocationScale','nu',nu);
% KS-Test with t-Distribution
[~,p] = kstest(XGivenT2Rand,'CDF',XDist);
disp(p)

%% Part 4
pVals = zeros(1,1000);

for iKSTest = 1:1000
    % Draw 10,000 samples from T2 Dist
    T2Rand = T2Dist.random(1,50);
    % For the Given T2 sample, draw corresponding X|T2 sample
    for iSample = 1:50
        XGivenT2Dist = makedist('Normal','mu',0,'sigma',sqrt(1/T2Rand(iSample))); 
        XGivenT2Rand(iSample) = XGivenT2Dist.random();
    end
    [~,pVals(iKSTest)] = kstest(XGivenT2Rand,'CDF',XDist);
end
% Plot Histogram of Values
pValHistFig = figure;
hist(gca,pVals,1000);
title('Histogram of p-values','Fontsize',14);
xlabel('Values','Fontsize',14);
ylabel('Count','Fontsize',14);
matlab2tikz(XHistFig,'pValHist.tikz');
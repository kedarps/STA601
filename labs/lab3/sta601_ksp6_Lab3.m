%% STA 601: Lab 3
% Author: Kedar S Prabhudesai
% Created on: 09/18/2013

close all;
clear all;

%% Full Conditionals
% Initial values
Beta0 = 0.05;

% N|beta,y ~ Poisson(25(1-beta))
NGivenBetaAndY = makedist('Poisson','lambda',25*(1-Beta0));
N1 = NGivenBetaAndY.random() + 20;
% beta|N,y ~ Beta(21,N-19)
BetaGivenNAndY = makedist('Beta','a',21,'b',N1-19);
Beta1 = BetaGivenNAndY.random();
% Samples from Full Conditionals
nSamples = 10000;
NSamples = zeros(1,nSamples);
BetaSamples = zeros(1,nSamples);

NSamples(1) = N1;
BetaSamples(1) = Beta1;

for iSample = 2:nSamples
    NGivenBetaAndY.lambda = 25*(1-BetaSamples(iSample-1));
    NSamples(iSample) = NGivenBetaAndY.random()+20;
    
    BetaGivenNAndY.b = NSamples(iSample)-19;
    BetaSamples(iSample) = BetaGivenNAndY.random();
end
% figure;plot(BetaSamples,NSamples,'bo-','Linewidth',2);
% title('First 10 Samples from Gibbs Sampler','FontSize',14);
% xlabel('\beta','FontSize',14);
% ylabel('N','FontSize',14);

% Burn-In
BetaSamples(1:1000) = [];
NSamples(1:1000) = [];

PostCredIntrval = quantile(BetaSamples,[0.05 0.95]);

figure;plot(1:numel(BetaSamples),BetaSamples,'bo','Linewidth',2);hold on;
plot(1:numel(BetaSamples),repmat(mean(BetaSamples),1,numel(BetaSamples)),'m-','Linewidth',3);
plot(1:numel(BetaSamples),repmat(PostCredIntrval(1),1,numel(BetaSamples)),'r-','Linewidth',3);
plot(1:numel(BetaSamples),repmat(PostCredIntrval(2),1,numel(BetaSamples)),'r-','Linewidth',3);hold off;
title('Trace Plots of \beta Samples post burn-in (1000)','FontSize',14);
xlabel('Iteration Number','FontSize',14);
ylabel('\beta','FontSize',14);

% P(N=20)
ProbOfInterest = mean(NSamples==20);

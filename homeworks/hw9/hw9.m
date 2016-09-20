%% STA 601 - Homework 9
% Author: Kedar Prabhudesai
% Created on: 10/09/2013

close all;
clear all;

%% This is the hieracrchical Model
% Yij ~ Exp(Li)
%  Li ~ Gamma(a,b)
%   a ~ Gamma(Aa,Ba)
%   b ~ Gamma(Ab,Bb)

%% Simulate Data
% Dummmy Variables to simulate data
nConditions = 10;
nSubjectsPerCondition = [100 75 65 100 80 90 100 95 70 60];
LambdaI = [0.3 0.1 1.5 0.5 0.7 2 0.9 1.1 0.2 2.1];
Yij = cell(1,nConditions);

% Create Distribution Objects
for iCond = 1:nConditions
    DistObj = makedist('Gamma','a',1,'b',1/LambdaI(iCond));
    Yij{1,iCond} = DistObj.random(1,nSubjectsPerCondition(iCond));
end

%% Gibbs Sampler
nGibbs = 5000;
nBurnIn = 1000;
LiSamples = zeros(nConditions,nGibbs);
aSamples = zeros(1,nGibbs);
bSamples = zeros(1,nGibbs);
%   p(a) ~ Exp(-a*Aa)
Aa = 5;
%   b ~ Gamma(Ab,Bb)
Ab = 4;
Bb = 5;
% Initialize
aSamples(1) = 2;
bSamples(1) = 3;
LiGivenAll = makedist('Gamma','a',2,'b',1/3);
bGivenAll = makedist('Gamma','a',Ab,'b',1/Bb);
SumYj = zeros(1,nConditions);

for iCond = 1:nConditions
    SumYj(iCond) = sum(Yij{1,iCond});
end

aGrid = (1:0.1:60);
aGridWeights = -aGrid.*Aa;

for iGibbs = 2:nGibbs
    home;
    disp(iGibbs);
    % Update Li for all i
    for iCond = 1:nConditions
        LiGivenAll.a = aSamples(iGibbs-1) + nSubjectsPerCondition(iCond);
        LiGivenAll.b = 1/(bSamples(iGibbs-1) + SumYj(iCond));
        LiSamples(iCond,iGibbs) = LiGivenAll.random();
    end
    
    % Update a 
    aProbWeights = (aGrid-1)*sum(log(LiSamples(:,iGibbs))) + aGridWeights;
    aProbWeights = exp(aProbWeights - max(aProbWeights));
    aSamples(iGibbs) = randsample(aGrid,1,true,aProbWeights);
    
    % Update b
    bGivenAll.b = 1/(Bb + sum(LiSamples(:,iGibbs)));
    bSamples(iGibbs) = bGivenAll.random();
end

% Burn-In
aSamples(1:nBurnIn) = [];
bSamples(1:nBurnIn) =  [];
LiSamples(:,1:nBurnIn) = [];

%% Manage Plotting
MeanLi = zeros(nConditions,1);
ConfInts = zeros(nConditions,2);
for iCond = 1:nConditions
    MeanLi(iCond) = mean(LiSamples(iCond,:));
    ConfInts(iCond,:) = quantile(LiSamples(iCond,:),[0.025 0.975]);
end

figure;
errorbar(1:nConditions,MeanLi,ConfInts(:,1),ConfInts(:,2),'Marker','diamond','Linewidth',2);hold on;
title('Point and Interval Estimates for different conditions','FontSize',14);
xlabel('Experiment Condition','FontSize',14);
ylabel('Mean Reaction Rate','FontSize',14);
ylim([0 8]);

plot(1:nConditions,LambdaI,'m*','Linewidth',2);hold off;
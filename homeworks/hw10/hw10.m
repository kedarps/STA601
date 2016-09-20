%% STA 601 - Homework 10
% Author: Kedar Prabhudesai
% Created on: 10/18/2013

close all;
clear all;

% Get data
tmp = importdata('china500.dat');
data = tmp.data;
data(:,[1 6]) = [];

% Outcome variable
y = data(:,1);
% Predictors
x = data(:,2:end);
x = [ones(size(x,1),1) x];
p = size(x,2);


%% Metropolis Hastings Algorithm
nSamples = 10000;
delta2 = 10;
betaSamples = zeros(p,nSamples);
% Predictor Weights. Initialize with random values
beta = mvnrnd(zeros(p,1),delta2*eye(p));
beta = beta';
betaSamples(:,1) = beta;

for iSample = 2:nSamples
    betaStar = mvnrnd(betaSamples(:,iSample-1),delta2*eye(p));
    betaStar = betaStar';
    
    home;
    disp(iSample);
    for iPred = 1:p
        betaZeroHere = betaSamples(:,iSample-1);
        betaStarHere = betaZeroHere;
        betaStarHere(iPred) = betaStar(iPred);
        
        LStarHere = sum(y.*log((1./(1+exp(-x*betaStarHere))) + (1-y).*log((1./(1+exp(x*betaStarHere))))));
        LZeroHere = sum(y.*log((1./(1+exp(-x*betaZeroHere))) + (1-y).*log((1./(1+exp(x*betaZeroHere))))));
        
        piStarHere = log(mvnpdf(betaStarHere,betaZeroHere,delta2*eye(p)));
        piZeroHere = log(mvnpdf(betaZeroHere,betaZeroHere,delta2*eye(p)));
        
        lhr = exp(LStarHere + piStarHere - (LZeroHere + piZeroHere));
        AcceptRejectFlag = binornd(1,min(1,lhr));
        
        if AcceptRejectFlag
            betaSamples(iPred,iSample) = betaStar(iPred);
        else
            betaSamples(iPred,iSample) = betaSamples(iPred,iSample-1);
        end
    end
end

% Manage Plotting
figure;
MeanBi = zeros(p,1);
ConfInts = zeros(p,2);
for iPred = 1:p
    MeanBi(iPred) = mean(betaSamples(iPred,:));
    ConfInts(iPred,:) = quantile(betaSamples(iPred,:),[0.025 0.975]);
end
figure;
errorbar(1:p,MeanBi,ConfInts(:,1),ConfInts(:,2),'Marker','diamond','Linewidth',2);
title('Point and Interval Estimates for Beta','FontSize',14);
xlabel('Predictor Index','FontSize',14);
ylabel('Value','FontSize',14);


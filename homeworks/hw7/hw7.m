%% STA 601 - Homework 7
% Author: Kedar Prabhudesai
% Created on: 9/24/2013

close all;
clear all;

%% Setup Data and Distributions
% Create Bivariate Distribution 
nSamples = 100;
rho = 0.8;
mu = [0 0]; 
SIGMA = [1 rho; rho 1];
y1 = -3:0.1:3;y2 = -3:0.1:3;
[Y1,Y2] = meshgrid(y1,y2);
yPDF = mvnpdf([Y1(:) Y2(:)],mu,SIGMA);
yPDF = reshape(yPDF,length(Y2),length(Y1));

% Get Data
rng('shuffle');
rSamples = mvnrnd(mu,SIGMA,nSamples);

% Find MLE
muMLE = mean(rSamples);
sigmaMLE = cov(rSamples);

%% Gibbs Sampler
nGibbs = 10000;
ySamples = rSamples';
ybar = mean(rSamples)';
mu0 = [0.2 0.2]';
L0 = [1.25 0.6;0.6 1.25];

nu0 = 4;
S0 = [1.2 0.4;0.4 1.3];

thetaSamples = zeros(2,nGibbs);
sigmaSamples = zeros(2,2,nGibbs);
sigmaSamples(:,:,1) = S0;

for iSample = 2:nGibbs
    % Update theta
    Ln = inv(inv(L0) + nSamples.*inv(sigmaSamples(:,:,iSample-1)));
    mun = Ln*(inv(L0)*mu0 + nSamples.*inv(sigmaSamples(:,:,iSample-1))*ybar);
    thetaSamples(:,iSample) = mvnrnd(mun,Ln);
    
    % Update Sigma
    Sn = S0 + (bsxfun(@minus,ySamples,thetaSamples(:,iSample)))*(bsxfun(@minus,ySamples,thetaSamples(:,iSample)))';
    Z = mvnrnd([0 0],inv(Sn),nu0+nSamples);
    sigmaSamples(:,:,iSample) = inv(Z'*Z);
end

% Burn-In
thetaSamples(:,1:1000) = [];
sigmaSamples(:,:,1:1000) = [];

% Estimate mu and Sigma
muBayes = mean(thetaSamples,2)';
sigmaBayes = mean(sigmaSamples,3);

%% Prediction
nPred = 50;
rng('shuffle');

muPred = mu;
sigmaPred = SIGMA;

yTruth = mvnrnd(mu,SIGMA,nPred);

y1PredTruth = zeros(nPred,1);
y1PredMLE = zeros(nPred,1);
y1PredBayes = zeros(nPred,1);

% Estimate (y1|y2) using truth
for iPred = 1:nPred
    y1PredTruth(iPred,1) = normrnd(muPred(1) + (sigmaPred(1,1)*sigmaPred(1,2)*(yTruth(iPred,2)-muPred(2)))/sigmaPred(2,2),sqrt((1-sigmaPred(1,2)^2)*sigmaPred(1,1)^2));
end
Ey1Giveny2Truth = mean(y1PredTruth)

% Estimate (y1|y2) using MLE
for iPred = 1:nPred
    y1PredMLE(iPred,1) = normrnd(muMLE(1) + (sigmaMLE(1,1)*sigmaMLE(1,2)*(yTruth(iPred,2)-muMLE(2)))/sigmaMLE(2,2),sqrt((1-sigmaMLE(1,2)^2)*sigmaMLE(1,1)^2));
end
Ey1Giveny2MLE = mean(y1PredMLE)

% Estimate (y1|y2) using Bayes Estimates
for iPred = 1:nPred
    y1PredBayes(iPred,1) = normrnd(muBayes(1) + (sigmaBayes(1,1)*sigmaBayes(1,2)*(yTruth(iPred,2)-muBayes(2)))/sigmaBayes(2,2),sqrt((1-sigmaBayes(1,2)^2)*sigmaBayes(1,1)^2));
end
Ey1Giveny2Bayes = mean(y1PredBayes);

MSPredErrMLE = (Ey1Giveny2Truth - Ey1Giveny2MLE).^2
MSPredErrBayes = (Ey1Giveny2Truth - Ey1Giveny2Bayes).^2

%% MLE Comparison
rng('shuffle');
yTruthMLEJoint = mvnrnd(muMLE,sigmaMLE,nPred);
EyTruthMLEJoint = mean(yTruthMLEJoint(:,1))

%% Confidence intervals for different predictors
ConfIntervalsTruth = quantile(y1PredTruth,[0.025 0.975])
ConfIntervalsMLE = quantile(y1PredMLE,[0.025 0.975])
ConfIntervalsBayes = quantile(y1PredBayes,[0.025 0.975])


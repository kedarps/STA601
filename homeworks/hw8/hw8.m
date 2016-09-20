%% STA 601 - Homework 8
% Author: Kedar Prabhudesai
% Created on: 9/27/2013

close all;
clear all;

%% Setup Data and Distributions
% Create Bivariate Distribution
nSamples = 100;
rho = 0.8;
mu = [0 0];
SIGMA = [1 rho; rho 1];

% Get Data
rng('shuffle');
yTruth = mvnrnd(mu,SIGMA,nSamples);

pctMissing = 0:0.05:0.5;

y1Bayes = [];
y1BayesConf = [];
y2Bayes = [];
y2BayesConf = [];

%% Missing Data
for iPct = 1:numel(pctMissing)
    % Simulate Missing Data
    yTruthMissing = yTruth;
    yTruthMissing = yTruthMissing';
    
    % Get Random indices to reject data
    OInds = rand(size(yTruthMissing));
    % Indicator which tells us what is missing: 1-Miss, 0-Present
    Oij = (OInds <= pctMissing(iPct));
    disp(['Processing ',num2str(pctMissing(iPct)),' Actual Percent ',num2str(sum(Oij(:))/numel(OInds))]);
    
    % Reject Data
    yTruthMissing(Oij) = NaN;
    % Make sure we do not have both y1 and y2 missing
    y1AndY2Missing = isnan(yTruthMissing(1,:)) & isnan(yTruthMissing(2,:));
    
    for iSample = 1:nSamples
        idxRnd = randperm(2,1);
        yTruthMissing(idxRnd,iSample) = yTruth(iSample,idxRnd);
        Oij(idxRnd,iSample) = 0;
    end
    idxMissing = find(sum(Oij,1)==1);
    
    %% Gibbs Sampler
    nGibbs = 10000;
    nBurnIn = 1000;
    
    mu0 = [0.2 0.2]';
    L0 = [1.25 0.6;0.6 1.25];
    
    nu0 = 4;
    S0 = [1.2 0.4;0.4 1.3];
    
    thetaSamples = zeros(2,nGibbs);
    sigmaSamples = zeros(2,2,nGibbs);
    
    % Initialize first value of Sigma and missing data
    sigmaSamples(:,:,1) = S0;
    for iMiss = 1:numel(idxMissing)
        % b: Missing Index
        b = find(Oij(:,idxMissing(iMiss)));
        yTruthMissing(b,idxMissing(iMiss)) = 2;
    end
    
    for iGibbs = 2:nGibbs
        % Update ybar
        ybar = mean(yTruthMissing,2);
        
        % Update theta
        Ln = inv(inv(L0) + nSamples.*inv(sigmaSamples(:,:,iGibbs-1)));
        mun = Ln*(inv(L0)*mu0 + nSamples.*inv(sigmaSamples(:,:,iGibbs-1))*ybar);
        thetaHere = mvnrnd(mun,Ln);
        thetaSamples(:,iGibbs) = thetaHere;
        % Update Sigma
        Sn = S0 + (bsxfun(@minus,yTruthMissing,thetaSamples(:,iGibbs)))*(bsxfun(@minus,yTruthMissing,thetaSamples(:,iGibbs)))';
        Z = mvnrnd([0 0],inv(Sn),nu0+nSamples);
        sigmaHere = inv(Z'*Z);
        sigmaSamples(:,:,iGibbs) = sigmaHere;
        
        % Update Missing Data
        for iMiss = 1:numel(idxMissing)
            % a: Present Index
            % b: Missing Index
            b = find(Oij(:,idxMissing(iMiss)));
            if b == 1
                a = 2;
            else
                a = 1;
            end
            ThetaMissGivenPres = thetaHere(b) + (sigmaHere(b,a)*(yTruthMissing(a,idxMissing(iMiss))-thetaHere(a)))/sigmaHere(a,a);
            SigmaMissGivenPres = sqrt(sigmaHere(b,b) - (sigmaHere(b,a)*sigmaHere(a,b))/sigmaHere(a,a));
            yTruthMissing(b,idxMissing(iMiss)) = ThetaMissGivenPres + SigmaMissGivenPres*randn;
        end
    end
    
    % Burn-In
    thetaSamples(:,1:nBurnIn) = [];
    sigmaSamples(:,:,1:nBurnIn) = [];
    
    % Estimate mu and Sigma
    muBayes = mean(thetaSamples,2)';
    sigmaBayes = mean(sigmaSamples,3);
    
    % Organize to plot
    y1Bayes = [y1Bayes;muBayes(1)];
    y1BayesConf = [y1BayesConf;quantile(thetaSamples(1,:),[0.025 0.975])];
    y2Bayes = [y2Bayes;muBayes(2)];
    y2BayesConf = [y2BayesConf;quantile(thetaSamples(2,:),[0.025 0.975])];
end
figure;
errorbar(pctMissing,y1Bayes,y1BayesConf(:,1),y1BayesConf(:,2),'Marker','diamond','Linewidth',2);
title('E(\theta_1) with 95% Credible Intervals','FontSize',14);
xlabel('p% Data Missig','FontSize',14);
% ylim([min(y1BayesConf(:,1))-0.1 max(y1BayesConf(:,2))+0.1]);
ylim([-1 1]);
xlim([-0.1 0.6]);

figure;
errorbar(pctMissing,y2Bayes,y2BayesConf(:,1),y2BayesConf(:,2),'Marker','diamond','Linewidth',2);
title('E(\theta_2) with 95% Credible Intervals','FontSize',14);
xlabel('p% Data Missig','FontSize',14);
% ylim([min(y2BayesConf(:,1))-0.1 max(y2BayesConf(:,2))+0.1]);
ylim([-1 1]);
xlim([-0.1 0.6]);

%% Complete Case Analysis
y1Bayes = [];
y1BayesConf = [];
y2Bayes = [];
y2BayesConf = [];

for iPct = 1:numel(pctMissing)
    % Simulate Missing Data
    yTruthMissing = yTruth;
    yTruthMissing = yTruthMissing';
    
    % Get Random indices to reject data
    OInds = rand(size(yTruthMissing));
    % Indicator which tells us what is missing: 1-Miss, 0-Present
    Oij = (OInds <= pctMissing(iPct));
    disp(['Processing ',num2str(pctMissing(iPct)),' Actual Percent ',num2str(sum(Oij(:))/numel(OInds))]);
    
    % Reject Data
    yTruthMissing(Oij) = NaN;
    % Make sure we do not have both y1 and y2 missing
    y1AndY2Missing = isnan(yTruthMissing(1,:)) & isnan(yTruthMissing(2,:));
    
    for iSample = 1:nSamples
        idxRnd = randperm(2,1);
        yTruthMissing(idxRnd,iSample) = yTruth(iSample,idxRnd);
        Oij(idxRnd,iSample) = 0;
    end
    idxMissing = find(sum(Oij,1)==1);
    
    %% Get rid of samples where data is missing
    for iMiss = 1:numel(idxMissing)
        yTruthMissing(:,idxMissing(iMiss)) = [0;0];
    end
    
    %% Gibbs Sampler
    nGibbs = 10000;
    nBurnIn = 1000;
    
    mu0 = [0.2 0.2]';
    L0 = [1.25 0.6;0.6 1.25];
    
    nu0 = 4;
    S0 = [1.2 0.4;0.4 1.3];
    
    thetaSamples = zeros(2,nGibbs);
    sigmaSamples = zeros(2,2,nGibbs);
    
    % Initialize first value of Sigma
    sigmaSamples(:,:,1) = S0;
    ybar = mean(yTruthMissing,2);
    
    for iGibbs = 2:nGibbs
        % Update theta
        Ln = inv(inv(L0) + nSamples.*inv(sigmaSamples(:,:,iGibbs-1)));
        mun = Ln*(inv(L0)*mu0 + nSamples.*inv(sigmaSamples(:,:,iGibbs-1))*ybar);
        thetaHere = mvnrnd(mun,Ln);
        thetaSamples(:,iGibbs) = thetaHere;
        
        % Update Sigma
        Sn = S0 + (bsxfun(@minus,yTruthMissing,thetaSamples(:,iGibbs)))*(bsxfun(@minus,yTruthMissing,thetaSamples(:,iGibbs)))';
        Z = mvnrnd([0 0],inv(Sn),nu0+nSamples);
        sigmaHere = inv(Z'*Z);
        sigmaSamples(:,:,iGibbs) = sigmaHere;
    end
    
    % Burn-In
    thetaSamples(:,1:nBurnIn) = [];
    sigmaSamples(:,:,1:nBurnIn) = [];
    
    % Estimate mu and Sigma
    muBayes = mean(thetaSamples,2)';
    sigmaBayes = mean(sigmaSamples,3);
    
    % Organize to plot
    y1Bayes = [y1Bayes;muBayes(1)];
    y1BayesConf = [y1BayesConf;quantile(thetaSamples(1,:),[0.025 0.975])];
    y2Bayes = [y2Bayes;muBayes(2)];
    y2BayesConf = [y2BayesConf;quantile(thetaSamples(2,:),[0.025 0.975])];
end
figure;
errorbar(pctMissing,y1Bayes,y1BayesConf(:,1),y1BayesConf(:,2),'Marker','diamond','Linewidth',2);
title('Complete Case Analysis - E(\theta_1) with 95% Credible Intervals','FontSize',14);
xlabel('p% Data Missig','FontSize',14);
% ylim([min(y1BayesConf(:,1))-0.1 max(y1BayesConf(:,2))+0.1]);
ylim([-1 1]);
xlim([-0.1 0.6]);

figure;
errorbar(pctMissing,y2Bayes,y2BayesConf(:,1),y2BayesConf(:,2),'Marker','diamond','Linewidth',2);
title('Complete Case Analysis - E(\theta_2) with 95% Credible Intervals','FontSize',14);
xlabel('p% Data Missig','FontSize',14);
% ylim([min(y2BayesConf(:,1))-0.1 max(y2BayesConf(:,2))+0.1]);
ylim([-1 1]);
xlim([-0.1 0.6]);
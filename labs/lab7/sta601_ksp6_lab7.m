%% STA 601: Lab 7
% Author: Kedar S Prabhudesai
% Created on: 11/1/2013

close all;
clear all;

% Get Data
X = importdata('data.txt');
% Prior Parameters
mu0 = 0;
t0 = 0.05;
a = 1;
% Remember that 'b' parameter in matlab's gamma function is in fact '1/b'
b = 0.05;
n = numel(X);

% Target Distribution for mu
muFullCond =  @(t,mu) exp(-0.5*t*sum((log(X)-mu).^2) - 0.5*t0*(mu-mu0)^2);

% Number of Trials
nTrials = 500;
% Burn-In
nBurnIn = 100;
% Proposal Distribution Std Dev
SCand = 0.5;

% Initialize
muSamples = zeros(1,nTrials);
tSamples = zeros(1,nTrials);

% Full Conditional distribution for tau
tFullCond = makedist('Gamma','a',a+n/2,'b',1/(0.5*sum((log(X)-muSamples(1)).^2) + b));
tSamples(1) = tFullCond.random();

% Gibbs Sampling
for iTrial = 2:nTrials
    home;
    disp(iTrial);
    % Update mu | tau,x using M-H
    
    % Step 1: Sample from mu'|mu(s)
    muPrime = normrnd(muSamples(iTrial-1),SCand);
    
    % Step 2: Compute Acceptance Ratio
    r = muFullCond(tSamples(iTrial-1),muPrime)/muFullCond(tSamples(iTrial-1),muSamples(iTrial-1));
    
    % Step 3: Accept/Reject
    u = rand;
    if u < r
        muSamples(iTrial) = muPrime;
    else
        muSamples(iTrial) = muSamples(iTrial-1);
    end
    
    % Update tau | mu,x 
    tFullCond.b = 1/(0.5*sum((log(X)-muSamples(iTrial)).^2) + b);
    tSamples(iTrial) = tFullCond.random();
end

% Burn-In
muSamples(1:nBurnIn) = [];
tSamples(1:nBurnIn) = [];

% Convert to Sigma^2
s2Samples = 1./tSamples;

% Manage Plotting
figure('Position',[67   304   922   345]);
plot(muSamples,'b-');
xlabel('Iterations','FontSize',14);
ylabel('Samples','FontSize',14);
title('\mu Samples','FontSize',14);

figure('Position',[67   304   922   345]);
plot(s2Samples,'b-');
xlabel('Iterations','FontSize',14);
ylabel('Samples','FontSize',14);
title('\sigma^2 Samples','FontSize',14);

% Find estimates of mean and variance
MeanFromSamples = exp(muSamples + s2Samples./2);
VarFromSamples = (exp(s2Samples) - 1).*exp(2.*muSamples + s2Samples);

EstMean = mean(MeanFromSamples);
MeanConfInts = quantile(MeanFromSamples,[0.025 0.975]);

EstVar = mean(VarFromSamples);
VarConfInts = quantile(VarFromSamples,[0.025 0.975]);

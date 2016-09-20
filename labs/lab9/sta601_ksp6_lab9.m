%% STA 601: Lab 8
% Author: Kedar S Prabhudesai
% Created on: 11/22/2013
function sta601_ksp6_lab9

close all;
clear all;

%% Get Data
X = importdata('data.txt');
X = X.data;

%% Prior Parameters
mu10 = 0.6;mu20 = 0.4;
t10 = 400;t20 = 400;
a1 = 1;a2 = 1;
% Remember that 'b' parameter in matlab's gamma function is in fact '1/b'
b1 = 0.05;b2 = 0.05;
n = numel(X);

% Target Distribution for Mu1 and Mu2
Mu1FullCond =  @(Xwd,t1,mu1) exp(-0.5*t1*sum((log(Xwd)-mu1).^2) - 0.5*t10*(mu1-mu10)^2);
Mu2FullCond =  @(Xwe,t2,mu2) exp(-0.5*t2*sum((log(Xwe)-mu2).^2) - 0.5*t20*(mu2-mu20)^2);

%% Gibbs Sampler
nGibbs = 12000;
nBurnIn = 2000;
% Proposal Distribution Std Dev
SCand = 0.5;

% Initialize
ThetaSamples = zeros(1,nGibbs);
ThetaSamples(1) = 5/7;
Mu1Samples = zeros(1,nGibbs);
Mu2Samples = zeros(1,nGibbs);

Tau1Samples = zeros(1,nGibbs);
Tau2Samples = zeros(1,nGibbs);
Tau1Samples(1) = 0.05;
Tau2Samples(1) = 0.05;

kVals = zeros(1,nGibbs);

for iGibbs = 2:nGibbs
    home;disp(iGibbs);
    % Draw W from Binomial(n,Theta)
    W = rand(n,1);
    W = W < ThetaSamples(iGibbs-1);
    
    % Split data according to weekdays/weekends
    XWkDays = X(W == 1);
    XWkEnds = X(W == 0);
    
    kVals(iGibbs-1) = sum(W);
    %% Update Theta | ---
    ThetaSamples(iGibbs) = betarnd(5+kVals(iGibbs-1),2+n-kVals(iGibbs-1));
    
    %% Update Mu1 | --- using M-H
    % Step 1: Sample from Mu1'|Mu1(s)
    Mu1Prime = normrnd(Mu1Samples(iGibbs-1),SCand);
    
    % Step 2: Compute Acceptance Ratio
    r = Mu1FullCond(XWkDays,Tau1Samples(iGibbs-1),Mu1Prime)/Mu1FullCond(XWkDays,Tau1Samples(iGibbs-1),Mu1Samples(iGibbs-1));
    
    % Step 3: Accept/Reject
    u = rand;
    if u < r
        Mu1Samples(iGibbs) = Mu1Prime;
    else
        Mu1Samples(iGibbs) = Mu1Samples(iGibbs-1);
    end
    
    %% Update Mu2 | --- using M-H
    % Step 1: Sample from Mu2'|Mu2(s)
    Mu2Prime = normrnd(Mu2Samples(iGibbs-1),SCand);
    
    % Step 2: Compute Acceptance Ratio
    r = Mu2FullCond(XWkEnds,Tau2Samples(iGibbs-1),Mu2Prime)/Mu2FullCond(XWkEnds,Tau2Samples(iGibbs-1),Mu2Samples(iGibbs-1));
    
    % Step 3: Accept/Reject
    u = rand;
    if u < r
        Mu2Samples(iGibbs) = Mu2Prime;
    else
        Mu2Samples(iGibbs) = Mu2Samples(iGibbs-1);
    end
    
    %% Update Tau1 | ---
    Tau1Samples(iGibbs) = gamrnd(a1+kVals(iGibbs-1)/2,1/(0.5*sum((log(XWkDays)-Mu1Samples(iGibbs)).^2)+b1));
    
    %% Update Tau2 | ---
    Tau2Samples(iGibbs) = gamrnd(a2+(n-kVals(iGibbs-1))/2,1/(0.5*sum((log(XWkEnds)-Mu2Samples(iGibbs)).^2)+b2));
end

%% Burn-In
ThetaSamples(1:nBurnIn) = [];
Mu1Samples(1:nBurnIn) = [];
Mu2Samples(1:nBurnIn) = [];
Tau1Samples(1:nBurnIn) = [];
Tau2Samples(1:nBurnIn) = [];
kVals(1:nBurnIn) = [];

% Convert Taus to sigma^2
S1Samples = 1./Tau1Samples;
S2Samples = 1./Tau2Samples;
keyboard
end
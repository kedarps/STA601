%% STA 601: Lab 8
% Author: Kedar S Prabhudesai
% Created on: 11/8/2013
function sta601_ksp6_lab8
close all;
clear all;

% Get Data
X = importdata('data.txt');
XWkDays = [];
XWkEnds = [];

% Separate data into Weekdays and Weekends
for iData = 2:numel(X)
    quotes = find(X{iData} == '"');
    day = X{iData}(quotes(1)+1:quotes(2)-1);
    value = str2double(X{iData}(1:quotes(1)-2));
    
    if strcmp(day,'Saturday') || strcmp(day,'Sunday')
        XWkEnds = cat(1,XWkEnds,value);
    else
        XWkDays = cat(1,XWkDays,value);
    end
end

[mu1Samples,s21Samples,xWkDSamples,EstMu1,Mu1ConfInts,EstS21,S21ConfInts] = GibbsSampler(XWkDays,0.6,400,1,0.05);

% Manage Plotting
figure('Position',[67   304   922   345]);
plot(mu1Samples,'b-');
xlabel('Iterations','FontSize',14);
ylabel('Samples','FontSize',14);
title('\mu_1 Samples','FontSize',14);

figure('Position',[67   304   922   345]);
plot(s21Samples,'b-');
xlabel('Iterations','FontSize',14);
ylabel('Samples','FontSize',14);
title('\sigma_1^2 Samples','FontSize',14);

[mu2Samples,s22Samples,xWkESamples,EstMu2,Mu2ConfInts,EstS22,S22ConfInts] = GibbsSampler(XWkEnds,0.2,400,1,0.05);

% Compute Probabilities of Interest
M1GtM2 = mean(mu1Samples > mu2Samples);
S1GtS2 = mean(sqrt(s21Samples) > sqrt(s22Samples));
X1GtX2 = mean(xWkDSamples > xWkESamples);

keyboard
    function [muSamples,s2Samples,xSamples,EstMu,MuConfInts,EstS2,S2ConfInts] = GibbsSampler(data,mu0,t0,a,b)
%         % Prior Parameters: mu0, t0, a, b
%         % Remember that 'b' parameter in matlab's gamma function is in fact '1/b'
        n = numel(data);
        
        % Target Distribution for mu
        muFullCond =  @(t,mu) exp(-0.5*t*sum((log(data)-mu).^2) - 0.5*t0*(mu-mu0)^2);
        
        % Number of Trials
        nTrials = 12000;
        % Burn-In
        nBurnIn = 2000;
        % Proposal Distribution Std Dev
        SCand = 0.5;
        
        % Initialize
        muSamples = zeros(1,nTrials);
        tSamples = zeros(1,nTrials);
        xSamples = zeros(1,nTrials);
        
        % Full Conditional distribution for tau
        tFullCond = makedist('Gamma','a',a+n/2,'b',1/(0.5*sum((log(data)-muSamples(1)).^2) + b));
        tSamples(1) = tFullCond.random();
        
        % Gibbs Sampling
        for iTrial = 2:nTrials
            home;disp(iTrial);
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
            tFullCond.b = 1/(0.5*sum((log(data)-muSamples(iTrial)).^2) + b);
            tSamples(iTrial) = tFullCond.random();
            
            % Get Samples from likelihood for Posterior Predictive
            xSamples(iTrial) = lognrnd(muSamples(iTrial),sqrt(1/tSamples(iTrial)));
        end
        
        % Burn-In
        muSamples(1:nBurnIn) = [];
        tSamples(1:nBurnIn) = [];
        xSamples(1:nBurnIn) = [];
        
        % Convert to Sigma^2
        s2Samples = 1./tSamples;
        
        % Get Estimates and Credible Intervals
        EstMu = mean(muSamples);
        MuConfInts = quantile(muSamples,[0.025 0.975]);
        
        EstS2 = mean(s2Samples);
        S2ConfInts = quantile(s2Samples,[0.025 0.975]);
    end

end
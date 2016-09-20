%% STA 601: Lab 6
% Author: Kedar S Prabhudesai
% Created on: 10/23/2013

close all;
clear all;

% Target Distribution
TarDist =  @(t) exp(-0.5.*(t.^2)) + 0.5*exp(-0.5.*((t-3).^2));
% Number of Trials
nTrials = 10000;
% Burn-In
nBurnIn = 2000;
% Proposal Distribution Std Dev
SCand = 100;
% Find number of Accepted Samples
nAccept = 0;

% Initialize
ThetaSamples = zeros(1,nTrials);

% Run Metropolis-Hastings
for iTrial = 2:nTrials
    home;
    disp(iTrial);
    % Step 1: Sample from theta'|theta(s)
    ThetaPrime = normrnd(ThetaSamples(iTrial-1),SCand);
    
    % Step 2: Compute Acceptance Ratio
    r = TarDist(ThetaPrime)/TarDist(ThetaSamples(iTrial-1));
    
    % Step 3: Accept/Reject
    u = rand;
    if u < r
        ThetaSamples(iTrial) = ThetaPrime;
        if iTrial >= nBurnIn
            nAccept = nAccept + 1;
        end
    else
        ThetaSamples(iTrial) = ThetaSamples(iTrial-1);
    end
end

% Compute Acceptance Ratio
AccRat = nAccept/numel(ThetaSamples);

% Theta Support to find analytic distribution
ThetaSupport = -4:0.01:8;
ThetaAnalytic = TarDist(ThetaSupport);

disp(AccRat);
figure;
axes = plotyy(ThetaSupport,ThetaAnalytic,ThetaSupport,ThetaAnalytic);
hold on
hist(axes(1),ThetaSamples,100);
ylim(axes(1),'Auto');
set(axes(1),'YTickMode','auto');
set(axes(2),'YTickMode','auto');
set(axes,'FontSize',14);

XLabelHandles = get(axes,'XLabel');
set(XLabelHandles{1},'String','\theta','FontSize',14);
YLabelHandles = get(axes,'YLabel');
set(YLabelHandles{1},'String','Count','FontSize',14);
set(YLabelHandles{2},'String','\pi(\theta)','FontSize',14);
LineHandle = get(axes(2),'Children');
set(LineHandle,'LineWidth',4);
title(['\sigma_{cand} = ',num2str(SCand)]);
hold off

%% STA 601 - Homework 3
% Author: Kedar Prabhudesai
% Created on: 9/11/2013

close all;
clear all;

%% Make Distributions
% Gamma Prior Parameters
a = 2; b = 1;
% Data for Women with Bachelors degree
n = 44; 
sumy = 66;
% Posterior
Posterior = makedist('Gamma','a',a+sumy,'b',1/(b+n));
% True mean of the Posterior
MeanFromPosterior = Posterior.mean();
% Probability of Interest from Distribution P(theta2 < 1.5)
ProbFromPosterior = Posterior.cdf(1.5);
% Lower and Upper 95% Bounds from Posterior
PostCredIntervals = Posterior.icdf([0.025 0.975]);

% Support for pdf
x = 0:0.01:4;
% Plot stuff
figure;
plot(x,Posterior.pdf(x),'Linewidth',2);hold on;
stem(PostCredIntervals,Posterior.pdf(PostCredIntervals),'Linewidth',2,'Marker','.','Color','b');
stem(MeanFromPosterior,0,'Linewidth',4,'Marker','*','Color','r');hold off;
title({'True Posterior',...
    ['Mean = ',num2str(MeanFromPosterior),' | 95% Credible Intervals = [',num2str(PostCredIntervals(1)),',',num2str(PostCredIntervals(2)),']'],...
    ['P(\theta_2 < 1.5) = ',num2str(ProbFromPosterior)]},...
    'FontSize',12);
xlabel('\theta | y','FontSize',12);
ylabel('p(\theta | y)','FontSize',12);
xlim([0.5 3.5]);

%% Do Monte-Carlo Simulations
nTrials = [10 100 1000];
ProbFromMC = zeros(size(nTrials));
MeanFromMC = zeros(size(nTrials));
CIFromMC = zeros(numel(nTrials),2);
rng('shuffle');

for iTrials = 1:numel(nTrials)
    % Draw Random Samples from Posterior
    PostSamples = Posterior.random(1,nTrials(iTrials));
    % Estimate the mean
    MeanFromMC = mean(PostSamples);
    % Estimate P(theta2 < 1.5)
    ProbFromMC = sum(PostSamples < 1.5)/nTrials(iTrials);
    % Find 95% Credible Intervals
    CredIntervalFromMC = quantile(PostSamples,[0.025 0.975]);
    
    % Plot stuff
    figure;
    plot(x,Posterior.pdf(x),'Linewidth',2);hold on;
    stem(CredIntervalFromMC,Posterior.pdf(CredIntervalFromMC),'Linewidth',2,'Marker','.');
    stem(MeanFromMC,0,'Linewidth',4,'Marker','*','Color','r');hold off;
    
    title({[num2str(nTrials(iTrials)),' MC Trials'],...
        ['Mean = ',num2str(MeanFromMC),' | 95% Credible Intervals = [',num2str(CredIntervalFromMC(1)),',',num2str(CredIntervalFromMC(2)),']'],...
        ['| MC Mean - True Mean | = ',num2str(abs(MeanFromPosterior-MeanFromMC))],...
        ['P(\theta_2 < 1.5) = ',num2str(ProbFromMC)]},...
        'FontSize',12);
    xlabel('\theta | y','FontSize',12);
    ylabel('p(\theta | y)','FontSize',12);
    xlim([0.5 3.5]);
end


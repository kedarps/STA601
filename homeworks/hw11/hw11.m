%% STA 601 - Homework 11
% Author: Kedar Prabhudesai
% Created on: 10/22/2013

close all;
clear all;

% Beta Prior Parameters
a = 1;
b = 1;
% Number of Trials
nTrials = 1:2000;
% We draw theta at random in [0,1]
H1Theta = rand(1,numel(nTrials));
% Get number of successes
k = binornd(nTrials,H1Theta);
% Calculate Bayes Factor
BayesFactor = beta(a+k,nTrials-k+b)./((0.5.^nTrials).*beta(a,b));

% Manage Plotting
figure;plot(nTrials,BayesFactor,'b-o');ylim([0 2]);hold on;
plot(nTrials,ones(numel(nTrials),1),'r','LineWidth',2);hold off;
xlabel('n','FontSize',14);
ylabel('\kappa','FontSize',14);
title(['\kappa vs Number of Trials, \alpha = ',num2str(a),', \beta = ',num2str(b)],'FontSize',14);
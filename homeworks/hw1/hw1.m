%% STA 601: Homework - 1

% Author: Kedar Prabhudesai
% Created on: 8/30/2013
close all;clear all;

% Number of students in class
n = 50;
% Expected Score - 82%
mu = 0.82;
% Std. Dev - 8%
sigma = 0.08;
% Support of the pdf
X = 0:0.01:1;
% Compute pdf over support
Y = (1/(sigma*sqrt(2*pi))).*exp(-((X-mu).^2)/(2*sigma^2));

% Make sure it is pdf. Verify that integration over support is 1
areaOverSupport = trapz(X,Y);
disp(['Area Under pdf = ',num2str(areaOverSupport)]);
% Plot Distribution
figure;plot(X,Y);
xlabel('Support (\theta)');
ylabel('p(\theta)');
title('Prior Density of \theta.');
%% Find 95% Confidence Interval

% Use Z* = 1.96
muUpper = mu + 1.96*sigma/sqrt(n);
muLower = mu - 1.96*sigma/sqrt(n);

disp(['95% Confidence Interval - ',num2str(muUpper),' < mu < ',num2str(muLower)]);


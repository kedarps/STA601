%% STA 601: Lab 5
% Author: Kedar S Prabhudesai
% Created on: 10/03/2013

close all;
clear all;

%% Gibbs Sampler with single variable update
% Initialize
mu1 = 0;
mu2 = [0;0];
S11 = 1;

nGibbs = 5000;
XSamples = zeros(1,nGibbs);
YSamples = zeros(1,nGibbs);
ZSamples = zeros(1,nGibbs);
% Starting values of Y, Z
YSamples(1) = 2;
ZSamples(1) = 3;

for iGibbs = 1:nGibbs-1
    % Update X
    S12 = [0.9 0.1];
    S21 = [0.9;0.1];
    S22 = [1 0.1;0.1 1];
    XSamples(iGibbs) = normrnd(mu1 + S12*inv(S22)*([YSamples(iGibbs);ZSamples(iGibbs)] - mu2),sqrt(S11 - S12*inv(S22)*S21));
    
    % Update Y
    YSamples(iGibbs+1) = normrnd(mu1 + S12*inv(S22)*([XSamples(iGibbs);ZSamples(iGibbs)] - mu2),sqrt(S11 - S12*inv(S22)*S21));
    
    % Update Z
    S12 = [0.1 0.1];
    S21 = [0.1;0.1];
    S22 = [1 0.9;0.9 1];
    ZSamples(iGibbs+1) = normrnd(mu1 + S12*inv(S22)*([XSamples(iGibbs);YSamples(iGibbs+1)] - mu2),sqrt(S11 - S12*inv(S22)*S21));
end

[XSamplesACF,lags] = autocorr(XSamples(1:1000),20);

figure('Position',[125 490 1175 400]);
subplot(1,2,1);
plot(XSamples(1:1000));
xlabel('Sample Number','Fontsize',14);
ylabel('X','Fontsize',14);
title('Sample Trace Plot','Fontsize',14);
subplot(1,2,2);
stem(lags,XSamplesACF);
xlim([-1 21]);
ylim([-1.1 1.1]);
xlabel('Number of Lags','Fontsize',14);
ylabel('Auto-Correlation','Fontsize',14);
title('Auto-Correlation Function','Fontsize',14);

%% Gibbs Sampler with Block Updates
XYSamples = zeros(2,nGibbs);
ZBlockSamples = zeros(1,nGibbs);
ZBlockSamples(1) = 5;

for iGibbs = 1:nGibbs-1
    % Update X,Y|Z
    mu1 = [0;0];
    mu2 = 0;
    S11 = [1 0.9;0.9 1];
    S12 = [0.1;0.1];
    S21 = [0.1 0.1];
    S22 = 1;
    XYSamples(:,iGibbs) = mvnrnd(mu1 + S12*inv(S22)*(ZBlockSamples(iGibbs) - mu2),sqrt(S11 - S12*inv(S22)*S21));
    
    % Update Z|X,Y
    mu1 = 0;
    mu2 = [0;0];
    S11 = 1;
    S12 = [0.1 0.1];
    S21 = [0.1;0.1];
    S22 = [1 0.9;0.9 1];
    ZBlockSamples(iGibbs+1) = normrnd(mu1 + S12*inv(S22)*(XYSamples(:,iGibbs) - mu2),sqrt(S11 - S12*inv(S22)*S21));
end

[XYSamplesACF,lags] = autocorr(XYSamples(1,1:1000),20);

figure('Position',[125 490 1175 400]);
subplot(1,2,1);
plot(XYSamples(1,1:1000));
xlabel('Sample Number','Fontsize',14);
ylabel('X','Fontsize',14);
title('Sample Trace Plot','Fontsize',14);
subplot(1,2,2);
stem(lags,XYSamplesACF);
xlim([-1 21]);
ylim([-1.1 1.1]);
xlabel('Number of Lags','Fontsize',14);
ylabel('Auto-Correlation','Fontsize',14);
title('Auto-Correlation Function','Fontsize',14);

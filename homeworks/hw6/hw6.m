%% STA 601 - Homework 6
% Author: Kedar Prabhudesai
% Created on: 9/19/2013

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

% Get MLE
muMLE = mean(rSamples);
sigmaMLE = cov(rSamples);

% Plot Contours
% figure;contour(y1,y2,yPDF,'Linewidth',2);
% colorbar;
% xlabel('y_1','Fontsize',14);
% ylabel('y_2','Fontsize',14);
% title('True Distribution','Fontsize',14);
% 
% figure;contour(y1,y2,reshape(mvnpdf([Y1(:) Y2(:)],muMLE,sigmaMLE),length(Y2),length(Y1)),'Linewidth',2);
% colorbar;
% xlabel('y_1','Fontsize',14);
% ylabel('y_2','Fontsize',14);
% title('Distribution from MLEs','Fontsize',14);
% 
% figure;contour(y1,y2,yPDF,'Linewidth',2);hold on;
% contour(y1,y2,reshape(mvnpdf([Y1(:) Y2(:)],muMLE,sigmaMLE),length(Y2),length(Y1)),'Linewidth',2);hold off;
% colorbar;
% xlabel('y_1','Fontsize',14);
% ylabel('y_2','Fontsize',14);
% title('Two Distributions Together','Fontsize',14);

%% Gibbs Sampler
nGibbs = 10000;
ySamples = rSamples';
ybar = muMLE';
mu0 = [0.2 0.2]';
L0 = [1.25 0.6;0.6 1.25];

nu0 = 4;
S0 = [1.2 0.4;0.4 1.3];
% S0 = [625 312.5;312.5 625];

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

muPost = mean(thetaSamples,2)';
sigmaPost = mean(sigmaSamples,3);
%% STA 601 - Homework 5
% Author: Kedar Prabhudesai
% Created on: 9/18/2013

close all;
clear all;

%% Initialize
% Total number of data 
n = 100;
% Number of Treated Subjects
m = 50;
% Draw sample data for y
yDist = makedist('Poisson','lambda',1);
ySamples = yDist.random(100,1);
% Also generate random values for x
xSamples = zeros(n,1);
xSamples(randperm(n,m)) = 1;

%% Using Gibbs Sampling
% Whether or not subject is treated
xs = 1;
% Total Number of Samples
nSamples = 10000;
% Initial value for Gamma
g0 = 4;
% Make Distribution Objects
LGivenGAndY = makedist('Gamma','a',sum(ySamples)+1,'b',1/(sum(g0.^xSamples)+n));
l1 = LGivenGAndY.random();
GGivenLAndY = makedist('Gamma','a',sum(ySamples.*xSamples)+1,'b',1/(l1*m+1));
g1 = GGivenLAndY.random();
YGivenLAndG = makedist('Poisson','lambda',l1*((1/g1)^xs));
y1 = YGivenLAndG.random();

% Samples for Full Conditionals
lSamples = zeros(nSamples,1);
gSamples = zeros(nSamples,1);
yPred = zeros(nSamples,1);

lSamples(1) = l1;
gSamples(1) = g1;
yPred(1) = y1;

% Perform Gibbs Sampling
for iSample = 2:nSamples
    LGivenGAndY.b = 1/(sum(gSamples(iSample-1).^xSamples)+n);
    lSamples(iSample) = LGivenGAndY.random();
    
    GGivenLAndY.b = 1/(lSamples(iSample)*m+1);
    gSamples(iSample) = GGivenLAndY.random();
    
    YGivenLAndG.lambda = lSamples(iSample)*((1/gSamples(iSample))^xs);
    yPred(iSample) = YGivenLAndG.random();
end

%% Find convergence diagnostics
% Autocorrelation
[lACF,lLags] = autocorr(lSamples,30);
[gACF,gLags] = autocorr(gSamples,30);

figure;stem(lLags,lACF);xlabel('Lags','Fontsize',14);ylabel('Autocorrelation','Fontsize',14);title('\lambda Autocorrelation function','Fontsize',14,'Linewidth',2);ylim([-1.2 1.2]);xlim([-1 31]);
figure;stem(gLags,gACF);xlabel('Lags','Fontsize',14);ylabel('Autocorrelation','Fontsize',14);title('\gamma Autocorrelation function','Fontsize',14,'Linewidth',2);ylim([-1.2 1.2]);xlim([-1 31]);

%% Estimate log(lambda) from Samples from Full Conditionals
logGammaEstimate = mean(log(gSamples));
Pct95CredInterval = quantile(log(gSamples),[0.025 0.975]);

disp(['Estimate: ',num2str(logGammaEstimate)]);
disp(['95 % Credible Intervals: ',num2str(Pct95CredInterval)]);
% Make Trace Plots
figure;plot(lSamples,gSamples,'bo');xlabel('\lambda','Fontsize',14);ylabel('\gamma','Fontsize',14);title('\lambda and \gamma Samples','Fontsize',14);
figure;plot(1:nSamples,gSamples,'bo');xlabel('Iteration Number','Fontsize',14);ylabel('\gamma','Fontsize',14);title('Samples from \gamma | \lambda,y','Fontsize',14);
figure;plot(1:nSamples,lSamples,'bo');xlabel('Iteration Number','Fontsize',14);ylabel('\lambda','Fontsize',14);title('Samples from \lambda | \gamma,y','Fontsize',14);

figure;plot(1:nSamples,log(gSamples),'bo');hold on;
plot(1:nSamples,repmat(logGammaEstimate,1,nSamples),'m-','Linewidth',2);
plot(1:nSamples,repmat(Pct95CredInterval(1),1,nSamples),'r-','Linewidth',2);
plot(1:nSamples,repmat(Pct95CredInterval(2),1,nSamples),'r-','Linewidth',2);hold off;
xlabel('Iteration Number','Fontsize',14);ylabel('log(\gamma)','Fontsize',14);title('Log of Samples from \gamma | \lambda,y','Fontsize',14);

figure;plot(1:nSamples,yPred,'bo');
xlabel('Iteration Number','Fontsize',14);ylabel('Y','Fontsize',14);title(['Samples from Predictive Posterior  Y|\lambda,\gamma for x^{(s)}=',num2str(xs)],'Fontsize',14);
figure;hist(yPred,10);
xlabel('Y','Fontsize',14);ylabel('Count','Fontsize',14);title(['Histogram of Samples from Predictive Posterior  Y|\lambda,\gamma for x^{(s)}=',num2str(xs)'],'Fontsize',14);
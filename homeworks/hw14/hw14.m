%% STA 601 - Homework 14
% Author: Kedar Prabhudesai
% Created on: 11/10/2013

close all;
clear all;

% Simulate Data
TrueBeta = [2 5];
% Predictors from Normal distribution with mean 1 and std dev. 2
X = 1 + 2.*randn(100,1);
% We append ones to X
X = cat(2,ones(100,1),X);
% Generate Z - Latent data
Z = X*TrueBeta' + randn(100,1);
% Set Y based on Z
Y = Z > 0;

XXInv = pinv(X'*X);

% Initialize prior values for bivariate beta prior
b0 = [0 0];
Sb = [0 0;0 0];

nTrials = 5000;
nBurnIn = 1000;
betaSamples = zeros(nTrials,2);
zSamples = zeros(nTrials,size(X,1));
% Initialize Latent data
zSamples(1,:) = rand(100,1);
zDistObj = makedist('Normal');

for iTrial = 2:nTrials
    home;disp(iTrial)
    % Update Beta
    bStar = XXInv*X'*zSamples(iTrial-1,:)';
    bHere = mvnrnd(bStar,XXInv);
    betaSamples(iTrial,:) = bHere;
    
    % Update z
    for iData = 1:size(X,1)
        xHere = X(iData,:);
        XB = xHere*bHere';
        
        zDistObj.mu = XB;
        zRand = zDistObj.random();
        
        if Y(iData) == 1
            while zRand < 0 
                zRand = zDistObj.random();
            end
        else
            while zRand > 0
                zRand = zDistObj.random();
            end
        end
        zSamples(iTrial,iData) = zRand;
    end
end

% Burn-In
betaSamples(1:nBurnIn,:) = [];
zSamples(1:nBurnIn,:) = [];
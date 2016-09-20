%% STA 601: Lab 1
% Author: Kedar S Prabhudesai
% Created on: 09/03/2013

% Support of pdf
x = 0:0.01:1;
%% Part 1: Find Posterior Distributions and Plot them
% Prior parameters
a = 0.5; b = 0.5;
% Likelihood Parameters
yA = 11; nA = 16;
yB = 5; nB = 6;

% Beta Distribution Object
pAPost = makedist('Beta','a',(a+yA),'b',(b+nA-yA));
pBPost = makedist('Beta','a',(a+yB),'b',(b+nB-yB));

% Plot Posterior
plot(x,pAPost.pdf(x),'b');hold on;
plot(x,pBPost.pdf(x),'r');hold off;
legend('Posterior for A','Posterior for B');

%% Part 2: 
% Find Probability of success at least greater than 0.8
% That is find area under the cdf greater than 0.8
pA = pAPost.icdf(1 - 0.8);
pB = pBPost.icdf(1 - 0.8);
disp(['p(A at least 0.8) = ',num2str(pA),'. p(B at least 0.8) = ',num2str(pB)]);

%% Part 3:
nTrials = 10000;
pARndVals = pAPost.random(1,nTrials);
pBRndVals = pBPost.random(1,nTrials);
TrueSuccessRate = sum(pBRndVals > pARndVals)/nTrials;

disp(['True Soln.B > Soln.A = ',num2str(TrueSuccessRate)]);
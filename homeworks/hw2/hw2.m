%% STA 601: Homework - 2
% Author: Kedar Prabhudesai; Created on: 09/04/2013

close all;clear all;

%% Part 1:
% We will use all the values of 'a' and 'b' and find how many trials we
% need in order to achieve Pr(theta <= 0.0015) = 0.95. We can find this
% probability by finding the value of cdf at theta = 0.0015.

% List of 'a' values of Prior
a = [1, 0.05, 1.6, 1.05];
% List of 'b' values of Prior
b = [666, 33.33, 407.4, 497];
% We are given that 'No Adverse Reactions occur'
y = 0;
% True value of theta
theta = 0.0015;

% Iterate through different Beta Priors
for iPrior = 1:numel(a)
    % Initial value of Number of Trials
    nTrials = 100;
    % Intial value of Pr(theta <= 0.0015)
    ProbOfInterest = 0;
    
    fprintf('Beta Prior Parameters, a = %f, b = %f\n',a(iPrior),b(iPrior));
    
    % We Increment Number of Trials until we achieve 
    % Pr(theta <= 0.0015) = 0.95
    while ProbOfInterest <= 0.95
        nTrials = nTrials + 1;
        % Calculate Posterior Beta parameters and make distribution
        aPosterior = (a(iPrior) + y);
        bPosterior = (b(iPrior) + nTrials - y);
        postDist = makedist('beta','a',aPosterior,'b',bPosterior);
        
        % Calculate the P(theta < 0.0015)
        % i.e. Calculate the value of cdf at 0.0015
        ProbOfInterest = postDist.cdf(theta);
    end
    fprintf('\tProbability = %f, Number Of Trials = %d\n\n',ProbOfInterest,nTrials);
end

%% Part 2:
% Now we will use a restricted prior Uniform[0 0.1]. 
% We can do this by choosing a Beta(1,1) prior, finding the posterior and
% then truncating it between (0,0.1)

% Initial value of Number of Trials
nTrials = 100;
% Intial value of Pr(theta <= 0.0015)
ProbOfInterest = 0;
% Uniform Prior Parameters
a = 1; b = 1;

fprintf('Beta Prior Parameters, a = %d, b = %d\n',a,b);

% We Increment Number of Trials until we achieve
% Pr(theta <= 0.0015) = 0.95
while ProbOfInterest <= 0.95
    nTrials = nTrials + 1;
    % Posterior Distribution is Binomial(n,p)
    postDist = makedist('Beta','a',(a + y),'b',(b + nTrials - y));
    % We truncate the posterior between (0,0.1)
    postDist = postDist.truncate(0,0.1);
    % Calculate the P(theta < 0.0015)
    % i.e. Calculate the value of cdf at 0.0015
    ProbOfInterest = postDist.cdf(theta);
end
fprintf('\tProbability = %f, Number Of Trials = %d\n\n',ProbOfInterest,nTrials);

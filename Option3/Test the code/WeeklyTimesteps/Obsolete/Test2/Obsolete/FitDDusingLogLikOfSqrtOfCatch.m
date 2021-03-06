% CREATED  18 April 2012
% MODIFIED 20 August 2013

% PURPOSE Fit a delay difference model to simulated data using
% a likelihood approach assuming the square root of catch in distributed as a Gaussian around the (square root) of predicted catch

% REFERENCE Application of a weekly delay-difference model to commercial catch and effort ...
%           Dichmont, Punt, Deng, Dell and Venables

% STATUS working

% BACKGROUND the delay difference model struggle to time series of CPUE if the starting biomasses Biomass(1) and (2) which were fixed are used to fit the time series of data. If nothing is done, the model will return un-realistic estimates to try to reconcile the model with the observation. The solution was to leave the first 2 years of data out of the objective function.

% Fitting the model against CPUEs to estimate catchability, 2 parameters for the distribution of recruitment and 10 total annual recruitment
% EKPdatasets; EKPparameters; [mle,fval,exitflag] = fminsearch(@(par) EKPLogLikOfCatch(par, cpue, true), [ones(1,13)],  optimset('MaxFunEvals',1e5, 'MaxIter', 1e5, 'Display','iter')); ConfidenceInterval2
 
% Plot observed CPUE against fitted:
%plot(cpue, '--rs'); hold; plot(105:520, 5*1e-5 * effort(105:end) ./ (5*1e-5 * effort(105:end) + M) .* Biomass(105:end) .* (1 - exp(- 5*1e-5 * effort(105:end) - M)) ./ effort(105:end)); hold off;

function negLL = EKPLogLikOfSqrtOfCatch(par, cpue, summation)

% ARGUMENT par is a vector of parameters
%          cpue is a time series of catch per unit of effort
%          sum is a boolean determining whether the function returns the sum of the negative log-likelihood or a vector of probability associated with each observation
global Biomass timesteps M ctch effort Recruitment RecSurvey Survival sigma

%%%%% Allocate individual parameters in vector par to clearer names
catchability_q = 5 * 1e-5;
% Compute the biomass dynamic
DelayDifference(par);

%% Calculate predicted cpue using only observation from the 105 month onward

% Using Quinn and Deriso (1999) equation for catch (Eq. 1.26 if I remember correctly)
pred_catch = catchability_q * effort(105:end) ./ (catchability_q * effort(105:end) + M) .* Biomass(105:end) .* (1 - exp(- catchability_q * effort(105:end) - M));

% Corresponding catches
catch_butfirst105 = ctch(105:timesteps);

% Assuming square root of catch are distributed as a Gaussian with mean square root of catch 
diff = log(catch_butfirst105) - log(pred_catch);

ss       = sum(diff.^2);
sss      = diff.^2;
nss = timesteps - 104;
sigma    = sqrt(ss./nss); % Estimate of the variance

cpueLL = nss * log(sqrt(2 * pi) * sigma) + 1/(2 * sigma.^2) * ss;
cpueLL_elements = nss * log(sqrt(2 * pi) * sigma) + 1/(2 * sigma.^2) .* sss;

% Same as above, expressed using Haddon (2001) simplification, Eq. 3.19 
% cpueLL = (nss./2).*(log(2*pi)+(2.*log(sigma))+1);% negative loglikelihood for a logNormal distribution 

if(summation == 1) negLL = cpueLL; end;
if(summation == 0) negLL = cpueLL_elements; end;

end


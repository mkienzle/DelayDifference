% CREATED  12 June 2012
% MODIFIED 12 June 2012

% PURPOSE calculate uncertainties on parameter estimates using the Jacobian as suggested by John D'Errico

% METHOD uses John D'Errico library called Adaptive Robust Numerical Differentiation available from MATLAB Central


% Create a function suitable for John D'Errico library
fun = @(par) EKPLogLikDelayDifference(par,cpue,0);

% Calculate the Jacobian at the Maximum Likelihood
[jac, err] = jacobianest(fun, mle);

% Calculate the covariance matrix (S) given the error variance (sigma) provided by EKPLogLikDelayDifference
global sigma
S = sigma * inv(jac' * jac);

% Parameter uncertainty
error = sqrt(diag(sigma * inv(jac' * jac)));
% CREATED  4 June 2012
% MODIFIED 4 June 2012

% PURPOSE estimate parameter uncertainty using the shape of the likelihood function and a certain (chi-square) distance from the minimum


%%%%% Plot the likelihood profile

% Number of points
n = 1e3;

% Results
results = [nan(n,1), nan(n,1)];

% Range of catchability
catchability = 0:20/n:20;

for i = 1:n

  results(i,1) = catchability(i);
  results(i,2) = EKPLogLikDelayDifference([catchability(i),x(2:13)], cpue);

end

% 50% Quantile of the chi-square distribution

Chisquare_quantile = 95.3;

hold on;
plot(results(:,1), results(:,2));
line([0; 20], [EKPLogLikDelayDifference([x], cpue) + Chisquare_quantile; EKPLogLikDelayDifference([x], cpue) + Chisquare_quantile]);
hold off;

%%%%% Search for the confidence limits

% Lower boundary
results2 = [nan(n,1), nan(n,1)];

% Range of catchability
tmp = 0:x(1)/n:x(1);

for i = 1:n

  results2(i,1) = tmp(i);
  results2(i,2) = (EKPLogLikDelayDifference([tmp(i),x(2:13)], cpue) - EKPLogLikDelayDifference([x], cpue) - Chisquare_quantile).^2;

end
[Y,I] = min(results2(:,2));
display(results2(I,1))

%plot(results2(:,1), results2(:,2))

% Upper boundary
results3 = [nan(n,1), nan(n,1)];

% Range of catchability
tmp2 = x(1) : (10*x(1) - x(1))/n : 10*x(1);

for i = 1:n

  results3(i,1) = tmp2(i);
  results3(i,2) = (EKPLogLikDelayDifference([tmp2(i),x(2:13)], cpue) - EKPLogLikDelayDifference([x], cpue) - Chisquare_quantile).^2;

end
[Y2,I2] = min(results3(:,2));
display(results3(I2,1))

%plot(results3(:,1), results3(:,2))
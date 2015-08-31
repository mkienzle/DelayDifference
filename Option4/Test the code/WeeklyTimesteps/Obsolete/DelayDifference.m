% CREATED  18 April 2012
% MODIFIED 20 Aug   2013

% PURPOSE create biomass dynamic according to the delay difference equation
% provided in Hilborn and Walters (1992) p.334, Eq. 9.2.11

function DelayDifference(par)

global Biomass Survival Recruitment timesteps effort M rho w_12weeksOld w_13weeksOld Tot_yr_rec sp_frac;

%%%%% Allocate individual parameters and scale them

% Catchability
%catchability_q = par(1) * 1e-5;
catchability_q = 5 * 1e-5;

% von Mises distribution parameters that distributes total annual recruitment to each month
%rbar = par(2);
%k = par(3);

% Total yearly recruitment for each years
%Tot_yr_rec(1:10,1) = par(2:11) * 1e9;
% NOTE NOTE NOTTE If you fix recruitment to value used in the simulation YOU GET A BIASED RESULT !
Tot_yr_rec(1:10,1) = [3.0 2.3 1.6 0.9 0.2 0.2 0.4 0.6 0.8 1.0] .* 1e9;
% Calculate recruitment pattern for the average fishing year
%rr=[1:12];%month(1:12); % where 1=Nov, ..., 12=Oct.
%deg=((rr./12)*(2*pi))'; %convert rr to radians of a circle
%rp=(2*pi*besseli(0,k))\exp(k.*cos(deg-(ones(12,1)*((rbar./12).*2.*pi))));%von mises distribution (Mardia and Jupp 2000 page 36) using matrix division operator \ rather than inv.
%phi_t=rp./(ones(12,1)*sum(rp)); 
%phi_t = transpose([0.3 0.15 0.1 0 0 0 0 0 0 0.05 0.1 0.3]);
phi_t = transpose([zeros(1,12) 1 zeros(1,39)]);
% Calculate survival using effort and natural mortality (M)
Survival = exp(-M) * exp( - catchability_q * effort);

% Distribute total recruitment in the first year over the first 12 months
Recruitment(1 : 52, 1) = Tot_yr_rec(1, 1) * phi_t;

%%%%% Dynamic of the biomass
for i=3:timesteps
    
    % Distribute recruitment over the next 12 monhts at the start of each year
    if(mod(i, 52) == 1)
        Recruitment(i : i+51, 1) = Tot_yr_rec(floor((i - 1) / 52) + 1, 1) * phi_t;
    end 
    
    % Calculate the biomass using the delay difference equation
    Biomass(i) = Survival(i-1) * Biomass(i-1) + rho * Survival(i-1) * Biomass(i-1) - rho * Survival(i-1) * Survival(i-2) * Biomass(i-2) - ...
        Survival(i-1) * rho * w_12weeksOld * Recruitment(i-1) + w_13weeksOld * Recruitment(i); 
    
end
end

% CREATED  18 April 2012
% MODIFIED 20 Aug   2013

% PURPOSE load data for delay difference model to be fitted
% to Eastern King Prawn

global ctch cpue effort
%%% Load simulated catch, effort and catch per unit of effort data (see SimulatePopDynamic.R)
% use only the last 10 years of data
ctch = csvread('Data/SimulatedCatch.csv', 1,1 );
ctch = ctch(2756 - 10 * 52 + 1 : 2756, 1);

%%% Load the last 10 years of cpue from a dataset containg monthly catch per boat days from 1958 to
%%% 2010 for all regions combined (a matrix 53 x 12)
cpue = csvread('Data/SimulatedCPUE.csv', 1,1 );
cpue = cpue(2756 - 10 * 52 + 1 : 2756, 1);

%%% Similarly load simulated effort
effort = csvread('Data/SimulatedEffort.csv', 1,1 );
effort = effort(2756 - 10 * 52 + 1 : 2756, 1);



 
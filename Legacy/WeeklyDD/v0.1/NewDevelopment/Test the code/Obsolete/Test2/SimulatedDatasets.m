% CREATED  18 April 2012
% MODIFIED 26 Aug   2013

% PURPOSE load data for delay difference model to be fitted
% to simulated data

global ctch effort
%%% Load simulated catch, effort and catch per unit of effort data (see SimulatePopDynamic.R)
% use only the last 10 years of data
ctch = csvread('Data/SimulatedCatch.csv', 1,1 );
ctch = ctch(2756 - 10 * 52 + 1 : 2756, 1);

%%% Similarly load simulated effort
effort = csvread('Data/SimulatedEffort.csv', 1,1 );
effort = effort(2756 - 10 * 52 + 1 : 2756, 1);

%%%%% Create a warm-up period of 10 years at the beginning of the dataset, made of a repetition of the first 52 weeks of data
A = repmat(ctch(1:52),10,1);
ctch = vertcat(A,ctch);

C = repmat(effort(1:52),10,1);
effort = vertcat(C,effort);


 
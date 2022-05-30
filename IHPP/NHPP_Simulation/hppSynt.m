function [tN] = hppSynt(n, lambda)
% Generate n exponentially-distributed
% inter-event times with mean mu

% interarrival times of Posson process are distributed 
% exponential with parameter lambda

% Simulating a Poisson process with rate Lambda
% (1) Generate ui?Uniform(0,1), i = 1, 2, ...,
    u = rand(n, 1);
% (2) Set yi = - ln(1-ui)/lambda
    y = -log(1-u)/lambda;
% (3) Set xi =sum(yj)
    tN = cumsum(y);
 end
function [ tK ] = nhppSynt( hLambda,params)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% T0 - time interval is defined [0, T0)
% L - is a function handle to rate funtion L(t)
% s - optional - set random generator seed
% Simulating a Poisson process with rate Lambda

% global params;

% if(exist >=3)
%     rng(s);
%     sprintf('Random generator seed set');
% end
global FIX_SEED
if FIX_SEED
    rng(params.seed);   
end
t=params.simTstart; 
tK=[];
x = linspace(params.simTstart,params.simTfinish,1000);
Lu = 10*max(hLambda(x)); %Fid HPP that setisfies Lu(t)>L(t) forevery t

while( t<params.simTfinish)   
    u = rand(1);   % (2) Generate ui~U(0,1), i = 1, 2, ...,    
    t = t - log(u)/Lu;    % (3) Next event if the PP(Lu)
    p = rand(1); % (4.1) Flip a coin
    if ( p<= hLambda(t)/Lu) % Event is not rejected
        tK=[tK; t];     %(4.2) Deliver event
    end
end

if (tK(end)>params.simTfinish)
    tK = tK(1:end-1);
end

end




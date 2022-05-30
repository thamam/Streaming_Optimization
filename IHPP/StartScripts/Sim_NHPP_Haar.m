%% Reconstructing Lamda(t) from sample of known arrival times
%% Single Batch Reconstruction
clear ;
% close all; 
clc
%% Initialize
% addpath('C:\Users\TOMER\Documents\Phd\Spring2016\ECE 6254\HW\HW6');
% s = rng;
%% Create fct. handle for L(t) - NHPP intensity function
% basisType = 'haar';
% hfunc = setBasisFunc(basisType);
J = 0;
T0 = 40;
Q =(T0*(1+J));
rate_mean = 3;
sigma = 2;
alpha = abs(rate_mean+sigma*randn(Q,1));
t=linspace(0,T0,1000);
% OW alpha with example values
% load('examples.mat');
% alpha = alpha_ex1;
% [f ] = baseFunc( t,alpha,J );
Type='piecewiseconst';
param = struct;
param.tr=[ 3  6  8 10 15 20 30];
param.val = [12 23 12 23 12 23 12 23 ];
lambda = generateIntensityFct(Type,T0,param,t);

%% Generate random sample point of L(t)
% L = @(t) baseFunc(t,alpha,J);
L = @(t) generateIntensityFct(Type,T0,param,t);
[ tN ] = nhppSynt( L,T0);
%% HPP - DEBUG
% rate=4;
% % [tN] = hppSynt(50,rate);
% % tN = tN(tN<T0);
% [tN] = HPPTRY(rate,T0);
% f= ones(size(t))*rate;

% Estimating Lambda(t) - i.e. intensity Function
Jhat = J;% = floor(log2(numel(tN)/T0));
d = ceil(T0/2^-Jhat); %number of basis in a single interval of [0,1] (resolution) times segmenth length in seconds
haarV0 = @(t) 1*((t>=0) & (t<1));
haarVj = @(t,n,j) 2^(j/2)*haarV0(2^j*t-n);
phi_n = @(t,n) haarVj(t,n,Jhat);
[alpha_hat,cntr,Gnorm] = negLogEstimate(phi_n,tN,T0,d);
cntr
%
tic
[alpha_hat_cvx] = lamdaEstCVX(phi_n,tN,T0,d);
toc
%
% Evaluating lambda_hat over all t in [0,T];
[Lambda ] = recLamda( t,alpha_hat,phi_n,d );
[Lambda_CVX ] = recLamda( t,alpha_hat_cvx,phi_n,d );

%
figure,
stairs(t,lambda,':g');title('\lambda(t)');xlabel('t');
hold all;
scatter(tN,L(tN),'*');
stairs(t,Lambda,'--r')
stairs(t,Lambda_CVX,'-+b')
ylabel(' \lambda^\^(t)');xlabel('t');
title(' Estimate of \lambda(t)');
legend('\lambda(t)','Events','\lambda^\^(t)','\lambda^\^_cvx(t)')
ylim([0,max(max(lambda),max(Lambda_CVX))*1.1]);
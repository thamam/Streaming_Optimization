%==========================================================================
%% Known bugs:
%==========================================================================
% 1. blksz==1 - prepare section fails to distribute basis function indices
% 2. argmincon with log barrier at the same time - the modified objective
% is used for both
% 3.TBD We need to select Time buffer and block size such that we can divide
% (tf-t0)/blocksize without a reminder
% 4.
%
%==========================================================================
%==========================================================================
%% Todo next:
%==========================================================================
% 1. Implement Equality constraints with Newton for the f_tild - i.e.
% change buffer evaluation mode when truncation is active Newton 2. 3.
%
%==========================================================================
%==========================================================================

%  Level 0
%% Title: Test bench for TONA
% Top simulation level that calls for the filtering algorithm
clear ;close all; clc;
%% Envrionment settings
addpath('KellyCODES/')
addpath(genpath('IHPP'))
% addpath('../NonLinearFiltering/')
addpath(genpath('auxFunc'));
addpath(genpath('Core'));
addpath(genpath('Utils'));
addpath(genpath('Optimizers'));

%% Primary Simulation Uset Input settings
bksz = 2; %block size
% TBD We need to select Time buffer and block size such that we can divide
% (tf-t0)/blocksize without a reminder
spord = 2;
MAXBUFSIZE = 10;

%Matlab optimization engine options
useObjGrad = false; % set true to use analytical gradient %%

%Set options for logbarrier function
UseLGBR = 0; %true; %set to '1' to use log barrier constraints
lb_delta = 1 ;  %weight starting value

% If nonegxMod==true, use inequality constraints to force x s.t. x>=0. If
% false, then we use the ineq constraints that x.'^psi(t)>=0 over a finite
% grid.
nonegxMod = 1;

% if runBatch==true then for each T, we compute the batch solution for
% [0,T] from scratch
runBatch=false;

% if runArgmincon==true, then argmincon solver will run and J_T W/O the log
% barrier part will be computed separatley.
runArgmincon = false;

streamOptions = struct('lb_delta',lb_delta,...
    'nonegxMod',nonegxMod,'bksz',bksz);

%% Simulation results and printout settings
% All printouts will be executed upon completion of simulation. Figures
% prints are detailed in PRINTOUTSCRIPTNAME script. Printing modes are :
%   1. PLOTFIGURES = 0;PLOTALL=*; Hold off all printing: 2. PLOTFIGURES =
%   1;PLOTALL=1; Print all posible figures. 3. PLOTFIGURES = 1;PLOTALL=0;
%   Print selectively. Set flags in PRINTOUTSCRIPTNAME
PRINTOUTSCRIPTNAME = 'main_nnlf_printouts.m';
PLOTFIGURES         = 0; % set to 0 to stop all printouts
PLOTALL             = 0; % set 1 to plot all possible figures (requires PLOTFIGURES==1)

%% Synthesizing the nonhomogenous Poisson process
% lambdaMin = 50;  lambdaMax = 200; t0 = 0 ; tf = 109;
lambdaMin = 5;  lambdaMax = 50; t0 = 0 ; tf = 43.5;
[tK, alphaGrndTr, nhpp_t , nhp_params ] = nhppsamp(lambdaMax, lambdaMin, spord, t0, tf );
fs = 100; % time resolution for discretization reconstruction
tvec = linspace(t0, tf, (tf-t0)*fs);
%% Run single TONA simulation
if 1
    [xstrm, XHIST, streambrkdwnData]  = TONA(nhp_params, streamOptions, tK, MAXBUFSIZE);
end
%% Run TONA W/O memory restriction
if 1
    [xstrm_unlim, XHIST_unlim, streambrkdwnData]  = TONA(nhp_params, streamOptions, tK, 2*numel(nhp_params.nInd));
end
%% Run simulation for different buffer sizes
if 1
    BUFSZARRAY = [1,2,3,4,5,6,7,8,9,10]*bksz;
    K = numel(BUFSZARRAY);
    ResArray = cell(numel(BUFSZARRAY),1);
    for k=1: K
        [ResArray{k}, XHIST_BUF]  = TONA(nhp_params, streamOptions, tK, BUFSZARRAY(k));
%         xrec  = xReconstruct(XHIST_BUF,bksz,BUFSZARRAY(k));
    end
    %   
end
%%
 
%%

%% Saving results from run with "infinite" memory
if false
    simParamsStr = sprintf('bksz_%i',bksz)
    resName = sprintf('SimResults_%s_%s.mat',date,simParamsStr);
    save(resName);
end


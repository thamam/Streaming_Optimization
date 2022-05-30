%==========================================================================
%% Top_Newton.m High level Summary:
%
%==========================================================================

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
% 2. Add highlevel summary
% 3. Add link in paper to code
% 4. Make a formal release branch on github
% 5. Add Justin to git hub
% 6. Remove from repository gim generated files (to save space)
% 7. Write readme file
% 8. Add more screen printouts
%==========================================================================

global DEBUG_MODE;
global FIX_SEED

DEBUG_MODE = false;
FIX_SEED = true;
%  Level 0
%% Title: Test bench for NOA
% Top simulation level that calls for the filtering algorithm
clear ;close all; clc;
% Envrionment settings
addpath('KellyCODES/')
addpath(genpath('IHPP'))
addpath(genpath('auxFunc'));
addpath(genpath('Core'));
addpath(genpath('Utils'));
addpath(genpath('Optimizers'));
addpath(genpath('PredCorr'));

%% Flags
flag_plotsplines = false;

%% Primary Simulation Uset Input settings
bksz = 1; %block size
% TBD We need to select Time buffer and block size such that we can divide
% (tf-t0)/blocksize without a reminder
spord = 2;
MAXBUFSIZE = 5;

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
lambdaMin = 50;  lambdaMax = 200; t0 = 0 ; tf = 109;
% lambdaMin = 5;  lambdaMax = 50; t0 = 0 ; tf = 43.5;
fprintf(' $$ Preparing simulation data \n\n');
[tK, alphaGrndTr, nhpp_t , nhp_params ] = nhppsamp(lambdaMax, lambdaMin, spord, t0, tf );
fs = 100; % time resolution for discretization reconstruction
tvec = linspace(t0, tf, (tf-t0)*fs);

if flag_plotsplines
    plot_basis_functions(nhp_params)
end
%% Run single NOA simulation
if 0
    fprintf(' $$ Handing over to Newton Online Algorithm (NOA) \n\n');
    [xstrm, XHIST]  = NOA(nhp_params, streamOptions, tK, MAXBUFSIZE);
    xtemprec = xReconstruct(XHIST, bksz, MAXBUFSIZE);
end
%% Prediction correction
if 1
    pc_bksz = 1;
    [xpc, xpc_iterhist]  = NOA_Prediction_Corrrection(nhp_params, streamOptions, tK, pc_bksz);

    if pc_bksz==1
        xpc_final = [diag(xpc,1);xpc(end)];
    else
        xpc_rec_mat = xReconstruct(xpc,pc_bksz,1);
        xpc_final = xpc_rec_mat(:,end);
    end
end
%% Run NOA W/O memory restriction
if 1
    [xstrm_unlim, XHIST_unlim, streambrkdwnData]  = NOA(nhp_params, streamOptions, tK, 2*numel(nhp_params.nInd));
end
%% Run simulation for different buffer sizes
if 0
    BUFSZARRAY = [1,2,3,4,5,6,7,8,9,10]*bksz;
    K = numel(BUFSZARRAY);
    ResArray = cell(numel(BUFSZARRAY),1);
    for k=1: K
        [ResArray{k}, XHIST_BUF]  = NOA(nhp_params, streamOptions, tK, BUFSZARRAY(k));
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



%%
% xpc_=[0;xpc];
lambda_pc = @(t) lambda_t(t, xpc_final ,nhp_params.phi_n, nhp_params.nInd);
lambda_nobuff = @(t) lambda_t(t,xstrm_unlim, nhp_params.phi_n, nhp_params.nInd);
lambda_NOA = @(t) lambda_t(t,xstrm, nhp_params.phi_n, nhp_params.nInd);

tvec = linspace(t0,tf,1e3);

figure(1),clf
hold all
plot(tvec,  lambda_NOA(tvec),'-.g');
% for i=2:2:size(XHIST,2)
%     Noa_t = @(t) lambda_t(t,XHIST(:,i), nhp_params.phi_n, nhp_params.nInd);
%     plot(tvec, Noa_t(tvec))
% end
plot(tvec, lambda_pc(tvec),'--r','linewidth',1.5)
plot(tvec,  lambda_nobuff(tvec),'k','LineWidth',1.5);

legend('NOA','pred-corr','NoBuff')
title('PC with bksz=1, solution takes after the second solve')
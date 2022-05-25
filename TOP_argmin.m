%==========================================================================
%% Known bugs:
%==========================================================================
% 1. blksz==1 - prepare section fails to distribute basis function indices
% 2. argmincon with log barrier at the same time - the modified objective
% is used for both 3.
% 
% 
% 
%==========================================================================
%==========================================================================

%==========================================================================
%% Todo next:
%==========================================================================
% 1. Implement Equality constraints with Newton for the f_tild - i.e.
% change buffer evaluation mode when truncation is active Newton 2. 3.
% 
% 
% 
%==========================================================================
%==========================================================================


%  Level 0
%% Example Title
% Top simulation level that calls for the filtering algorithm

%% Envrionment settings
clear all
close all
clc

addpath('KellyCODES/')
addpath(genpath('../NHPP/'))
addpath('../NonLinearFiltering/')
addpath(genpath('auxFunc'));
%

%% Primary Simulation Uset Input settings
bksz = 2; %block size
% We need to select Time buffer and block size such that we can divide
% (tf-t0)/blocksize without a reminder
spord = 2;
MAXBUFSIZE = 5;

%Matlab optimization engine options
%% Known bug - objective also includes logbarrier penalty
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

%% Synthesizing the nonhomogenous Poisson process
% lambdaMin = 50;  lambdaMax = 200; t0 = 0 ; tf = 103;
lambdaMin = 5;  lambdaMax = 50; t0 = 0 ; tf = 43.5;
[tK, alphaGrndTr, nhpp_t , nhp_params ] = nhppsamp(lambdaMax, lambdaMin, spord, t0, tf );


fs = 100; % time resolution for discretization reconstruction
tvec = linspace(t0, tf, (tf-t0)*fs);

%% Simulation results and printout settings
% All printouts will be executed upon completion of simulation. Figures
% prints are detailed in PRINTOUTSCRIPTNAME script. Printing modes are :
%   1. PLOTFIGURES = 0;PLOTALL=*; Hold off all printing: 2. PLOTFIGURES =
%   1;PLOTALL=1; Print all posible figures. 3. PLOTFIGURES = 1;PLOTALL=0;
%   Print selectively. Set flags in PRINTOUTSCRIPTNAME
PRINTOUTSCRIPTNAME = 'main_nnlf_printouts.m';
PLOTFIGURES         = 1; % set to 0 to stop all printouts
PLOTALL             = 0; % set 1 to plot all possible figures (requires PLOTFIGURES==1)

%% Initializing the streaming simulation settings
nhp_params.blocksize = bksz;
nhp_params.funcNum = nhp_params.d/bksz-1;
nhp_params.blocksNum = nhp_params.funcNum + 1;
nhp_params.fminconOptions = ...
    optimoptions(@fmincon,'MaxIterations',300,'MaxFunctionEvaluations',50000, ...
    'OptimalityTolerance',1e-10,'ConstraintTolerance',1e-10);

%additional solver options settings
nhp_params.fminconOptions.SpecifyObjectiveGradient = useObjGrad;


[Tlim_div, xidmapTbl , xIndTbl, xlimTbl] = streamingbreakdown(nhp_params);
Tlim_div=[0,0;Tlim_div]; % first row is zeros because T starts from 1 and blocks start from 0
Tend = size(Tlim_div,1);

%% Initialization
JTObj     = JT_BuffObjCLASS(MAXBUFSIZE, bksz);
% JTObj= JT_BuffObjCLASS(MAXBUFSIZE, bksz);
xbufObj   = MEMBUFFERCLASS(MAXBUFSIZE+1,bksz,'double',true); % buffer to manage active xvariables - \xtld_t
% size set to MAXBUFSIZE+1 to match function buffer which has a size of
% MAXBUFSIZE
JTinfObj  = cell(Tend,1);
% do not change data written to frist row of arracy - saved for info on
% array -
JTinfObj{1,2} = bksz;
x0=(1:1:bksz).'*0+1;
xhat_strm=x0;
Xflushed=[]; %xflushed records all the optimization variable that were discarded from the active memory(buffer)
Ts = (Tlim_div(end)-Tlim_div(2,1))/numel(tK);
% Ts = (Tlim_div(end)-Tlim_div(2,1))/30;
%% Main loop - simulating streaming arrival of batches
tic
for T=2:Tend
    JTinfObj{1,1} = T ; % 'Saved space - ptr to last row with with data
    
    %==========================================================================================%
    % Instantiate fT and add it to the array
    % %
    %==========================================================================================%
    xidmap = {xidmapTbl(T-1,:), xidmapTbl(T,:)};                                               %
    xInd ={xIndTbl(T-1,:), xIndTbl(T,:)};                                                      %
    
    nhpObj_T = NHP_ftObjCLASS(spord, tK, Tlim_div(T,:), xInd,lb_delta);   
    fT = ftObjCLASS(nhpObj_T,T,xidmap);    % Input - new increment of f_T(x_{T-1},x_T)         %               
    JTObj.addfunction(fT,T);    
    JTinfObj{T,1} = fT;
    JTinfObj{T,2} = T;
    %==========================================================================================%
    xT0=ones(fT.nt,1);
    if T>=MAXBUFSIZE+2
        x0=[xhat_strm(bksz+1:end);xT0];
    else
        x0=[xhat_strm;xT0];
    end        
    %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    %% Run streamlined solver
    options.nonegxMod = nonegxMod;

    options.Solver = 'NW'; 
    mu=10;
    eps=1e-6;

    [xhat_acon, fval_strm, y0, DBG_strm]   =  argmincon(JTObj,x0, ...
            Ts,nhp_params.fminconOptions, options);        
    %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    xhat_strm = xhat_acon;
    %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\
    
    %% Save xhat to hist array
    indShift = bksz*max(0,T-(MAXBUFSIZE));
    maxind = min((T),MAXBUFSIZE);
    histInd = indShift + (1:bksz*maxind);
    if T>MAXBUFSIZE
        XHIST_argmin(histInd,T) = xhat_strm(bksz+1:end);
    else
        XHIST_argmin(histInd,T) = xhat_strm;
    end
    
    
    %% Option for batch run - [0, Tcurrent]
    if runBatch
        %%
        [xbatch, fval_batch, DBG_batch] = batchSolver(Tlim_div,tK, xIndTbl,...
            T, xidmapTbl, nhpObj_T, spord, Ts, nhp_params, options);                
    end
    
    %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    if T>MAXBUFSIZE
        Xflushed=[Xflushed;xhat_strm(1:bksz)];
    end
    
end
toc
%% Computing refernce batch solution
if false
    [ alpha_hat_cvx,  Lambda_CVX ] = nhpp_opt_solver(tK, nhp_params, tvec, [] );
end

%% Show Results
if (false)
    %% Ptintout prepared for the UNBUFFERED case
    tvecExtended = linspace(t0, tf + nhp_params.splineOrder, (tf-t0)*fs);
    figure(1),clf
    tittext = sprintf('NHPP Simulation ');
    title(tittext)
    hold on
    plot(tvecExtended, nhpp_t(tvecExtended), ':b','LineWidth',2.0   );
    plot(tK,ones(numel(tK),1)*0.2,'+k');
    plot(tvecExtended, lambda_t( tvecExtended ,alpha_hat_cvx ,nhp_params.phi_n,nhp_params.nInd), '.g','LineWidth',2); %% if xMlSol is array need to change into xMlSol{end}
    plot(tvecExtended, lambda_t( tvecExtended ,xhat_strm ,nhp_params.phi_n,nhp_params.nInd), 'r','LineWidth',1)
    legend('Ground Truth','Events' , 'Global' , 'Filtering')
end

%% Ptintout prepared for the BUFFERED case

if (false)
    %%     xbuff = [Xhist;xhat_(bksz+1:end)];
    xbuff= xhat_strm;
    tvecExtended = linspace(t0, tf + nhp_params.splineOrder, (tf-t0)*fs);
    %     figure(1),clf
    figure
    tittext = sprintf('NHPP Simulation ');
    title(tittext)
    hold on
    plot(tvecExtended, nhpp_t(tvecExtended), ':b','LineWidth',2.0   );
    plot(tK,ones(numel(tK),1)*0.2,'+k');
    plot(tvecExtended, lambda_t( tvecExtended ,alpha_hat_cvx ,nhp_params.phi_n,nhp_params.nInd), '.g','LineWidth',2); %% if xMlSol is array need to change into xMlSol{end}
    plot(tvecExtended, lambda_t( tvecExtended ,xbuff ,nhp_params.phi_n,nhp_params.nInd), 'r','LineWidth',1)
    legend('Ground Truth','Events' , 'Global' , 'Filtering')
    
end



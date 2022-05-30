function [tK, alphaSyn, nhpp_t, params  ] = nhppsamp(lambdaMax, lambdaMin, spord, t0, tf  )
% nhppsamp Simulate nonhomodenous Poisson process observations and
% generating Bsplines parameters
% Input:lambdaMax - 
%       lambdaMin - 
%       spord - 
%       t0 - 
%       tf - 
% Output: tK - 
%         alphaSyn - 
%         nhpp_t - 
%         params: 
%                'Prate'-
%                'minRate' - 
%                'splineOrder' -
%                'simTstart' -
%                'simTfinish' -
%                'iterThs' -
%                'sectionLength' - 
%                'nInd' - 
%                'SupLen'
%                'tsup.Start' -
%                'tsup.End' - 
%                'd'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%
%% Set Global values
global DEBUG_MODE ;

%% Configure simulation Settings
params = struct('Prate',lambdaMax,'minRate', lambdaMin,'splineOrder',spord,'simTstart',t0, 'simTfinish',tf,...
    'iterThs',1e-4,'sectionLength',1);
% loading previous sim data options0
UseSavedNhppSimdData = true;

%% NHPP_PREPARE_DATA.m
% Initialize
[params  ] = NHPP_GCD_Initialize(params); %initializing function
fs = 1e4;
tt =  params.simTstart:1/fs:params.simTfinish;
params.tt=tt;
dateOfFile=date;

%% Synthesize NHPP simulation Data
% Try to use simulation samples from previous runs
    if UseSavedNhppSimdData
        try
            load(params.fileNhppName);
            fprintf('$ Simulation data file was loaded sucessfully \n');
        catch ME
            if( strcmp(ME.identifier,'MATLAB:load:couldNotReadFile') )
                ME.message
                fprintf('Simulating new process \n');
                [alphaSyn] = syntContIntNHPP(params); %generate coeffs for Lambda(t)
                [ tK ] = nhppSynt(@(t)lambda_t(t,alphaSyn,params.phi_n,params.nInd) ,params); %send function handle of lambda as function of t
                save(params.fileNhppName,'tK','alphaSyn','dateOfFile','params');
            end
        end
    else % Required new simulation
        fprintf('Simulating new process data \n');
        [alphaSyn] = syntContIntNHPP(params); %generate coeffs for Lambda(t)
        [ tK ] = nhppSynt(@(t)lambda_t(t,alphaSyn,params.phi_n,params.nInd) ,params); %send function handle of lambda as function of t
        if(OverWrite_NhppSimData)
            dateOfFile=date;
            save(params.fileNhppName,'tK','alphaSyn','dateOfFile');
        end
    end

%% Generate random sample point of L(t)
nhpp_t = @(t) lambda_t(t,alphaSyn,params.phi_n,params.nInd) ; 
%         MLE.theoretical = likelyHoodEstimateSpline(tK,params,alphaSyn);
%  [Lambda_GroundTruth]      = lambda_t(tt,alphaSyn,params.phi_n,params.nInd);
% [Lambda_at_arrivalPoints]   = lambda_t(tK,alphaSyn,params.phi_n,params.nInd);


end


%% Sim_NHPP_Bsplines
%% Reconstructing Lambda(t) from sample of known arrival times
%% Notation
% Lambda(t) - Ground truth intensity function
% CVX - sloves min log likelihood s.t. lambda(t)>=0
% DCM - solves by iterative method based on Nestrov's dual averaging
%%
clear ;close all; clc
addpath('DCM_Based_Approach','BasisFunctions','NHPP_Simulation',...
    'Coord_Descent');
%% Set Simulation params
RUNMode.CVX = true;
RUNMode.DCM = true;

lambdaMax = 5 ; %maximal rate for Lambda(t) 
global params
params = struct('Prate',lambdaMax,'splineOrder',3);
[params.simTstart, params.simTfinish ]  = deal (0,15);

UseSavedNhppSimdData = true;
OverWrite_NhppSimData = false;
Save_data = true;

opts.sectionLength = 5; opts.MaxIter =25;

%% Initialize
nhpp_Initialize(); %initializing function
fs = 1e2;
tt =  params.simTstart:1/fs:params.simTfinish;

if UseSavedNhppSimdData 
    try 
        load(params.fileNhppName);
    catch ME
        if( strcmp(ME.identifier,'MATLAB:load:couldNotReadFile') )
            ME.message
            fprintf('Simulating new process \n');
            [alphaSyn] = syntContIntNHPP(params); %generate coeffs for Lambda(t)
            [ tK ] = nhppSynt(@(t)lambda_t(t,alphaSyn,params.phi_n,params.nInd) ,params); %send function handle of lambda as function of t
            dateOfFile=date;
            save(params.fileNhppName,'tK','alphaSyn','dateOfFile');
        else 
            load(params.fileNhppName);
        end   
    end    
else % Required new simulation
    fprintf('Simulating new process \n');
    [alphaSyn] = syntContIntNHPP(params); %generate coeffs for Lambda(t)
    [ tK ] = nhppSynt(@(t)lambda_t(t,alphaSyn,params.phi_n,params.nInd) ,params); %send function handle of lambda as function of t
    if(OverWrite_NhppSimData)
        dateOfFile=date;
        save(params.fileNhppName,'tK','dateOfFile');
    end
end

%% Generate random sample point of L(t)
MLE.theoretical = likelyHoodEstimateSpline(tK,params,alphaSyn);
[Lambda_GroundTruth]        = lambda_t(tt,alphaSyn,params.phi_n,params.nInd);
[Lambda_at_arrivalPoints]   = lambda_t(tK,alphaSyn,params.phi_n,params.nInd);

% Plot theoretical lambda and it's sampled arrival times
figure(1),plot(tt,Lambda_GroundTruth,':g');title('\lambda(t)');xlabel('t'); hold all;
scatter(tK,Lambda_at_arrivalPoints,'*'); ylabel(' \lambda^\^(t)');xlabel('t'); xlabel('time');title(' Estimate of \lambda(t)'); legendVec={'\lambda(t)';'Events Samples'};
%% Run Coordinate Descent
tic
[alpha_hat_sec,sectionParams,tK_sec_ind,tSimSections]  = NHP_CoordDescent(tK,params,opts);
toc
%% Estimate Lambda(t) from samples using Dual averager solver
if RUNMode.DCM
    sprintf('starting DCM \n');
    [alpha_hat_dcm,iterNum] = dualAvgSolver(tK,params);   
    MLE.DCM         = likelyHoodEstimateSpline(tK,params,alpha_hat_dcm);
    [Lambda_DCM ]               = lambda_t(tt,alpha_hat_dcm,params.phi_n,params.nInd);    
    plot(tt,Lambda_DCM.','-.k'); legendVec=[legendVec;'DCM'];
end
%% Estimate Entire Sample using negative log likelihood minimization
if RUNMode.CVX
    [alpha_hat_cvx] = cvxEstLambdaSplines(tK,params);
    MLE.cvx = likelyHoodEstimateSpline(tK,params,alpha_hat_cvx);
    [Lambda_CVX ]               = lambda_t(tt,alpha_hat_cvx,params.phi_n,params.nInd);
    plot(tt,Lambda_CVX.','-b'); legendVec=[legendVec;'NLL Est'];
end
%% Plotting
legend(legendVec);

%%
saveRes = input('Would you like to save the results? (y/n)','s');
if strcmp(lower(saveRes),'y')
    save(['results_',date,params.fileNhppName]);
end

%% Printing Options

if false
    plotScript_coordniateDescent;
    
end



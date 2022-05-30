% SimUpdateStream

clear ;close all; clc
addpath('BasisFunctions','NHPP_Simulation',...
    'Coord_Descent','CVX_','BwdFwdDscnt','DisplayScripts','wavesProp');
%CVX_Support_Test();
% GLOBALS - UPPER CASE VARIABLES

global FIX_SEED
global DEBUG_MODE
% Set Global values
FIX_SEED = false;
DEBUG_MODE = false;

%% Set Sol Settings



%% Prepare Data Settings
lambdaMax = 25 ; %maximal rate for Lambda(t)
params = struct('Prate',lambdaMax,'splineOrder',3,'simTstart',2, 'simTfinish',59,...
    'iterThs',1e-4);
% loading previous sim data options0
UseSavedNhppSimdData = true;
OverWrite_NhppSimData = true;
Save_data = true;

%% Prepare Simulation Data
run('NHPP_PREPARE_DATA')
%%  Set Sovler Options
params.solverOptions = optimoptions(@fmincon,'MaxIterations',5000,'MaxFunctionEvaluations',50000,...
    'OptimalityTolerance',1e-5);

%%
for sim=SIM
    switch sim{:}.Type
        case 'GlobalCvx'
            try
                run('CVX_Global_Solution_SCRIPT');
            catch M
                Exception.CVX = M;
                fprintf('CVX Exception')
            end
        case 'CoordinateDescent'
            try
                %%
                sim{:}.options.alpha_hat_cvx = alpha_hat_cvx;
                [CoordDescOut]  = wavesCD(tK,params,sim{:}.options);
            catch M
                Exception.GCD = M;
                fprintf('GCD Exception\')
            end
        case 'wavesCD'
            try
                %
                sim{:}.options.alpha_hat_cvx = alpha_hat_cvx;
                [WavesCoordDescOut]  = wavesCD(tK,params,sim{:}.options);
            catch M
                Exception.WCD = M;
                fprintf('WCD Exception\n')
            end
        case 'FWD_SWEEP_W_TAIL'
            try
                %
                sim=SIM(4)
                sim{:}.options.alpha_hat_cvx = alpha_hat_cvx;
                [FWD__SWP__W__TAIL]  = wavesCD(tK,params,sim{:}.options);
            catch M
                Exception.BTBO = M;
                fprintf('BTBO Exception\n')
            end
        case 'SWEEP'
            try
                %
                sim=SIM(5);
                sim{:}.options.alpha_hat_cvx = alpha_hat_cvx;
                [SWEEPOut]  = wavesCD(tK,params,sim{:}.options);
            catch M
                Exception.SWEEP = M;
                fprintf('SWEEP Exception\n')
            end
            
        case 'RANDOMCD'
            try
                %
                sim=SIM(6);
                sim{:}.options.alpha_hat_cvx = alpha_hat_cvx;
                [RANDOMCDOut]  = wavesCD(tK,params,sim{:}.options);
            catch M
                Exception.RANDOMCD = M;
                fprintf('RANDOMCD Exception\n')
            end      
        otherwise
            error('Simulation Type is un recognizable \n');
    end
end
saveRes = 'y'; %input('Would you like to save the results? (y/n)','s');
if strcmpi(saveRes,'y')
    formatOut = 'D-mm_dd_yy_T-HH_MM_PM';
    timeStamp = datestr(datetime('now'),formatOut);
    save(['resultsSTRM_',timeStamp,'_',params.fileNhppName]);
end



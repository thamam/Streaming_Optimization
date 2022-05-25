%% SimUnified
%% Creadted on - June 13
%% SimulationParamsIteration
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
%% Configure simulation Settings
lambdaMax = 24 ; %maximal rate for Lambda(t)
params = struct('Prate',lambdaMax,'splineOrder',1,'simTstart',0, 'simTfinish',20,...
    'iterThs',1e-4,'sectionLength',1);
% loading previous sim data options0
UseSavedNhppSimdData = true;
OverWrite_NhppSimData = true;
Save_data = true;
%% Prepare Simulation Data
run('NHPP_PREPARE_DATA')
%% Set Simulation models to run
run('SIM_SOL_CELL_SETTEING')

%%  Set Sovler Options
params.solverOptions = optimoptions(@fmincon,'MaxIterations',300,'MaxFunctionEvaluations',50000,...
    'OptimalityTolerance',1e-5);
%% 


solList_obj = CellList([numel(SIM),1]);
for sim=SIM
    solList_obj.add({sim{:}.Type});
    if strcmpi(sim{:}.Type,'GlobalCvx')
        run('CVX_Global_Solution_SCRIPT');
        
        if 0 %save data for use in nonlinear filtering
            save('sim_Tsim_rate_25_0_12_sec_1_sp_1','alpha_hat_cvx','alphaSyn','Lambda_CVX','Lambda_GroundTruth','params','tK','Lambda_at_t')
        end
    else
        sim{:}.options.alpha_hat_cvx = alpha_hat_cvx;
        eval(sprintf('Out_%s = wavesCD(tK,params,sim{:}.options);', sim{:}.Type));
    end
end
saveRes = 'y'; %input('Would you like to save the results? (y/n)','s');
if strcmpi(saveRes,'y')
    formatOut = 'D-mm_dd_yy_T-HH_MM_PM';
    timeStamp = datestr(datetime('now'),formatOut);
    save(['resultsUNIF_',timeStamp,'_',params.fileNhppName]);
end
%%
% SectionIterationCount
% itrCntVec = cumsum(sum(Out_wavesCD.activeFramesMat(1:Out_wavesCD.SectionIterationCount,:),2))
nrm0 = norm( alpha_hat_cvx );
% diffMat = Out_wavesCD.alphaValM(1:Out_wavesCD.SectionIterationCount,:) - alpha_hat_cvx(:).';
% relresnrm = sqrt(sum(diffMat.^2,2)) / nrm0;

 %%
%  lambda_waves = lambda_t(tt,Out_wavesCD.alphaValM(Out_wavesCD.SectionIterationCount,:).',params.phi_n,params.nInd);
% norm(Lambda_CVX-Lambda_GroundTruth)/norm(Lambda_GroundTruth)
% norm(Lambda_CVX-lambda_waves)/norm(Lambda_CVX)

%% Run NHPP part
rho = 0.35;
rtol=1e-4;
P=2;
 [ cnzs_x, cnzs_it_hist, cnzs_outstat,dbgdata ]  = admmnhppsolver( rho,rtol,params,P,tK );

 %%
 diffzatk=0;
 for i=1:numel(cnzs_it_hist)
    diffzatk(i) = norm(dbgdata.z(:,i) - alpha_hat_cvx)/norm(alpha_hat_cvx);
 end
    lambda_admm = lambda_t(tt,cnzs_x,params.phi_n,params.nInd);
    diffzatk2 = sqrt(sum((dbgdata.z -alpha_hat_cvx).^2,1))/norm(alpha_hat_cvx,2);
%%
figure(1),clf
plot(tt,Lambda_CVX)
hold all
plot(tt,lambda_admm,'--g','linewidth',1.5)
legend('glbl','admm')
%%
itrvcntvec =1: Out_wavesCD.SectionIterationCount;
figure(2),clf
semilogy(itrvcntvec,relresnrm,'-+r')
hold all
semilogy(cnzs_outstat(:,1),diffzatk(:),'-*b')
legend('CD-Waves','ADMM')
xlabel('itr[k]')
ylabel('||x^^-X-glbl||/\\X_glbl||')
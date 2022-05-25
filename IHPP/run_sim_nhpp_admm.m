% %% ADMM_Consensus_form_NHPP
clear
% close all
clc
%% Load sim data
addpath(genpath('../../Code(SVN)/'))
load('GCD_simData_maxRate_25_T_2_59_L_1.mat');
% load GCD_simData_maxRate_20_T_0_21_L_2.mat;
% load('GCD_simData_maxRate_60_T_0_21_L_3.mat')
% load GCD_simData_maxRate_20_T_0_500_L_3;
%%
tic
nrmerr=[];
rtol=1e-4;
rho = 0.01;
%rho relaxes the constraints - rememeber that it is also the step size
%of the gradient escent of the dual variable
maxitr = 50;
parms(1) = maxitr;
parms(2)=rtol;
parms(3) = rho;
parms(4) = params.d; %equiavalent to N - size of global variable
%% Prepare admm consensus problem for least squares
% P=10;
params.rho = rho;
Pvec = [5,10,15,20,40,50,75,100];
% Pvec = [5,100];
results = cell(size(Pvec));
figure(2),
hold all
for i=1:numel(Pvec)
    P = Pvec(i);
    parms(5)=P; % # of agents
    [Fi,gmap,Zgind,InitValues,xupdateGrad] =  prepadmmNHPP(P,params,tK);
    %% Run ADMM iterative Updates
    [cnzs_x, cnzs_it_hist, cnzs_outstat,dbgdata]= ...
        nhpp_admmconsensus(Fi,InitValues,parms,gmap,Zgind,[],@ZupdaadmmNhpp,xupdateGrad);
    cnzs_errnrm = norm(cnzs_x - alphaSyn)/norm(alphaSyn) ;
    
   results{i}.cnzs_x=cnzs_x;
   results{i}.cnzs_it_hist=cnzs_it_hist;
   results{i}.dbgdata=dbgdata;
   results{i}.P=P;
   results{i}.cnzs_outstat=cnzs_outstat;
    
    
    %% Plot Residuals
    itrcnt = size(cnzs_outstat,1);
    plot(cnzs_outstat(2:end,1),cnzs_outstat(2:itrcnt,3),'-*');
    xlabel('iteration[k]')
    ylabel('$\mathbf{\frac{||z^{k+1}-z^{k}||}{||z^0||}}$','interpreter','Latex');
    title('Consensus relative residual norm  $\mathbf{\frac{||z^{k+1}-z^{k}||}{||z^0||}}$','interpreter','Latex');
    %%
end
    toc
%%
figure(3),clf
hold all
for i=1:numel(Pvec)   
    itrcnt = numel(results{i}.cnzs_it_hist);
    plot(1:itrcnt,results{i}.cnzs_it_hist/results{i}.cnzs_it_hist(1),'-*');
end
    xlabel('iteration[k]')
    ylabel('$\mathbf{\frac{||z^{k+1}-z^{k}||}{||z^0||}}$','interpreter','Latex');
    title('Consensus relative residual norm  $\mathbf{\frac{||z^{k+1}-z^{k}||}{||z^0||}}$','interpreter','Latex');
   legend(num2str(Pvec.'))
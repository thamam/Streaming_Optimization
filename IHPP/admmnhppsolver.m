function [ cnzs_x, cnzs_it_hist, cnzs_outstat,dbgdata ] = admmnhppsolver( rho,rtol,params,P,tK )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
nrmerr=[];
% rtol=1e-6;
%rho relaxes the constraints - rememeber that it is also the step size
%of the gradient escent of the dual variable
% rho = 0.01;
maxitr = 300;
parms(1) = maxitr;
parms(2)=rtol;
parms(3) = rho;
parms(4) = params.d; %equiavalent to N - size of global variable
%% Prepare admm consensus problem for least squares
% P=3;
parms(5)=P; % # of agents
params.rho = rho;
[Fi,gmap,Zgind,InitValues,xupdateGrad] =  prepadmmNHPP(P,params,tK);

%% Run ADMM iterative Updates
[cnzs_x, cnzs_it_hist, cnzs_outstat,dbgdata]= ...
    nhpp_admmconsensus(Fi,InitValues,parms,gmap,Zgind,[],@ZupdaadmmNhpp,xupdateGrad);
% cnzs_errnrm = norm(cnzs_x - alphaSyn)/norm(alphaSyn) ;

end


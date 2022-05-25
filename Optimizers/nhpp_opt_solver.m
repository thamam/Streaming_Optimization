function  [ alpha_hat_cvx,  Lambda_star ,WX] = nhpp_opt_solver(tK, params, tt, x0 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
WX=[];
%% CVX regression - minimize over entire set of samples

% [alpha_hat_cvx] =  cvxEstLambdaSplines(tK,params);
%%
[ alpha_hat_cvx ] =  cvxEstLambdaSplines2(tK,params, x0);
% CmpError = sqrt(sum(alpha_hat_cvx-alpha_hat_cvx2).^2)
%%
% MLE.cvx = likelyHoodEstimateSpline(tK,params,alpha_hat_cvx);
% Lambda_star is the value of the coeff assuming that all data-past and
% future is availalbe
[ Lambda_star ]    = lambda_t(tt , alpha_hat_cvx , params.phi_n , params.nInd );

Lambda_at_t = @(x) lambda_t(tt,x,params.phi_n,params.nInd);

T = params.funcNum;

if (false && T>10)
    for t=1:T
        WX_params.nInd= params.nlimDiv(t,1):params.nlimDiv(t+1,2);
        WX_params.d = 2*params.blocksize;
        WX_params.phi_n = params.phi_n;
        WX_params.solverOptions = params.solverOptions;
        WX_params.simTstart = params.Tlim_div(t,1);
        WX_params.simTfinish= params.Tlim_div(t,2);   
        WX_tK =   tK(tK>WX_params.simTstart & tK<WX_params.simTfinish ) ;
        WX_x0 = ones(WX_params.d,1);
        wstar(:,t) =   cvxEstLambdaSplines2(WX_tK,WX_params, WX_x0);
        
    end
    WX=wstar;
    
end

end


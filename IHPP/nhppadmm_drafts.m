
%% Run Direct MLE Optimization
if false
    params.solverOptions = optimoptions(@fmincon,'MaxIterations',50,'MaxFunctionEvaluations',50000,...
        'OptimalityTolerance',1e-3);
    [alpha_hat_cvx] =  cvxEstLambdaSplines2(tK,params);
end

if true
%% Compare Direct to ADMsM method
    %     relError = norm(alpha_hat_cvx - cnzs_x)/norm(alpha_hat_cvx); 
    %% Visualize
    if false
        tt= linspace(params.simTstart,params.simTfinish,1000);
        [Lambda_CVX ]               = lambda_t(tt,alpha_hat_cvx,params.phi_n,params.nInd);
        OptLambda  = lambda_t(tt,alphaSyn,params.phi_n,params.nInd);
        admmLambda =  lambda_t(tt,cnzs_x,params.phi_n,params.nInd);
        figure(1),clf
        hold all
        plot(tt,OptLambda)
        %     plot(tt,Lambda_CVX)
        plot(tt,admmLambda)
        legend('Optilmal','Global','ADMM')
        xlabel('t')
        ylabel('I(t)')
    end
   %% Error convergence Analysis
%    diffzatk = sqrt(sum((dbgdata.z -alpha_hat_cvx).^2,1))/norm(alpha_hat_cvx,2);
    kitr = size(cnzs_outstat,1);
    figure(2),
%     subplot(121),
    plot(2:kitr,cnzs_outstat(2:kitr,3),'-*r');
    xlabel('iteration[k]')
    ylabel('$\mathbf{\frac{||z^{k+1}-z^{k}||}{||z^0||}}$','interpreter','Latex');
    title('Consensus relative residual norm  $\mathbf{\frac{||z^{k+1}-z^{k}||}{||z^0||}}$','interpreter','Latex');
%     hold all;
%     plot(2:kitr,diffzatk(2:kitr),'-*b');
%     legend({'$z-res$','$zvs/z_{glb}$'},'Interpreter','latex') % - Notice the usage of {} to encompass all the legend entries together
    %     subplot(122),
%     semilogy(2:kitr,cnzs_outstat(2:kitr,3),'-+r');
%     hold all;
%     semilogy(2:kitr,diffzatk(2:kitr),'-*b');
%     xlabel('iteration[k]')
%     ylabel('$\mathbf{log \left( \frac{||z^{k+1}-z^{k}||}{||z^0||} \right)}$','interpreter','Latex');
%     title('Consensus relative residual norm  $\mathbf{log \left( \frac{||z^{k+1}-z^{k}||}{||z^0||} \right)}$','interpreter','Latex');
%     
%     legend({'$z-res$','$zvs/z_{glb}$'},'Interpreter','latex') % - Notice the usage of {} to encompass all the legend entries together
end

end

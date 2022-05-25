PLOTFIGURES =1;
prnt_DynSelUpdatesConvColumnSubfig =1;
prnt_UpdatesConvColumnSubfig=1;


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

%% Plot Updates conversion - Superimposed
% We show here the convergence of the updates by plotting the estiamte for
% four consecutive batches, starting from 'pind' . We then draw the
% attention of the reader to the estimate of the same x_t across all four
% estimates to see it convergence to the value of the batch estiamte - x*_t
if (true)
    
    lambda_hat = lambda_t( tvec ,alpha_ref ,nhp_params.phi_n,nhp_params.nInd);


    xsol = XHIST_unlim;
    legarr={};
    colArr = {'.b','-.g' ,'--r' ,'m'};
    pind =    13;
    lambda_estp0 = lambda_t( tvec ,xsol(:,pind)   ,nhp_params.phi_n,nhp_params.nInd);
    lambda_estp1 = lambda_t( tvec ,xsol(:,pind+1) ,nhp_params.phi_n,nhp_params.nInd);
    lambda_estp2 = lambda_t( tvec ,xsol(:,pind+2) ,nhp_params.phi_n,nhp_params.nInd);
    lambda_estp3 = lambda_t( tvec ,xsol(:,pind+3) ,nhp_params.phi_n,nhp_params.nInd);
    
    figure(5),clf
    tittext = sprintf('Visualizing updates conversion');
    title(tittext)
    hold all
    %plotting the refrence estimate
    plot(tvec, lambda_hat, 'k','LineWidth',2,'MarkerSize',8,'MarkerIndices',1:4:length(tvec));
    %plotting estimate at consecutive times
    plot(tvec, lambda_estp0, '.-','LineWidth',1.5,'MarkerSize',8,'MarkerIndices',1:4:length(tvec));
    plot(tvec, lambda_estp1, '.','LineWidth',1.5,'MarkerSize',8,'MarkerIndices',1:4:length(tvec));
    plot(tvec, lambda_estp2, '-','LineWidth',1.5,'MarkerSize',8,'MarkerIndices',1:4:length(tvec));
    plot(tvec, lambda_estp3, '--','LineWidth',1.5,'MarkerSize',8,'MarkerIndices',1:4:length(tvec));
    
    %
    legarr{1} = sprintf('$ \\lambda^*(t)$');
    
    for p=1:numel(pindArray)
        pind = pindArray(p);
        locnlim =sectionsPlan(pind,:) ;
        localnInd = params.nInd(params.nInd>=locnlim(1) & params.nInd<=locnlim(2));
        plot(tvec, lambda_t( tvec ,xsol(:,pind) ,params.phi_n,localnInd),colArr{p},'LineWidth',1.5,'MarkerIndices',1:3:length(tvec));
        legarr{p+1} = sprintf('$ \\hat{\\lambda_{%i}(t)}$', pind);
    end
    xlim([0 30])
    h_leg =  legend(legarr,'Interpreter','latex');
    h_leg.FontSize=12;
    h_leg.Location='northwest';
end
%         title( '$ \| x_{filt} - x_{glbl} \|   $  vs. Memory  Length ','Interpreter','latex')
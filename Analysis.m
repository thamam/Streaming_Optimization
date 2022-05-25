%% Envrionment settings
addpath('KellyCODES/')
addpath(genpath('../NHPP/'))
% addpath('../NonLinearFiltering/')
addpath(genpath('auxFunc'));

%% Plotting Mode flags
close all;
clc;
clear;

% set flag plotall=1 to plot all figures, or select selectively by changing
% the value in the if from 0 to 1 before each plot
plotall = false;

%% Load Results file
% load('SimResults_30-Jul-2020_bksz_2.mat');
load('SimResults_31-Jul-2020_bksz_2.mat');

%% Setting the batch solution as the reference
% [ alpha_hat_cvx,  Lambda_CVX ] = nhpp_opt_solver(tK, nhp_params, tvec, [] );
% alpha_ref = alpha_hat_cvx;
% load('xest_nomemlim.mat'); % will load 'XHIST_UNLIM'  variable
alpha_ref = XHIST_unlim(:,end);


%% Use this plot to compare single estimate from buffere reconstruction vs the batched solution
% XHIST shold be the results array from the buffered simulation
if (1 || plotall)
    %%
    xrec  = xReconstruct(XHIST,bksz,MAXBUFSIZE);
    % This gives the final solution with the fixed values retrieved into
    % xbuf_full
    xbuf_full = xrec(:,end);
    
    % The plot shows the batch estiamte reconstruction Vs. buffered
    % reconstruction
    figure(1),clf
    subplot(2,1,1)
    tvecExtended = linspace(t0, tf + nhp_params.splineOrder, (tf-t0)*fs);
    tittext = sprintf('NHPP Simulation ');
    title(tittext)
    hold on
    % plot(tvecExtended, nhpp_t(tvecExtended), ':b','LineWidth',2.0   );
    plot(tK,ones(numel(tK),1)*0.2,'+k');
    plot(tvecExtended, lambda_t( tvecExtended ,alpha_ref ,nhp_params.phi_n,nhp_params.nInd), '.g','LineWidth',2); %% if xMlSol is array need to change into xMlSol{end}
    plot(tvecExtended, lambda_t( tvecExtended ,xbuf_full ,nhp_params.phi_n,nhp_params.nInd), 'r','LineWidth',1)
    legend('Events' , 'Batch' , 'Buffered')
    
end
%% This figure plots the estimate error of \|x_t^* - x~_t\|_2 for all t<= T
if (0 || plotall)
    T = size(streambrkdwnData.Tlim_div,1);
    e_tld_2_variables = abs(xbuf_full - alpha_ref);
    e_tld_blocks = zeros(T,1);
    for t=1:(T+1)
        blkind = (t-1)*bksz + (1:bksz);
        xrec_t = xbuf_full(blkind);
        xhat_t = alpha_ref(blkind);
        e_tld_blocks(t) = norm(xrec_t - xhat_t,2);
    end
    %
    d = numel(e_tld_2_variables);
    tmp_rr = repmat(e_tld_blocks,2,1);
    blck_err_repmat = tmp_rr(:).';
    % figure(2), clf
    subplot(2,1,2)
    % subplot(211)
    hold all
    plot(1:d,e_tld_2_variables,'.-r')
    % subplot(212)
    plot(1:d,blck_err_repmat/bksz,'.-b') % plot of normalized block error
end


%% Plot Updates conversion - Column subfigures - Fixed time reference
% We show here the convergence of the updates by plotting the estiamte for
% four consecutive batches, starting from 'pind' . We then draw the
% attention of the reader to the estimate of the same x_t across all four
% estimates to see it convergence to the value of the batch estiamte - x*_t
if (1|| plotall)
    Tlim_div =    streambrkdwnData.Tlim_div;
    xbatchsol = XHIST_unlim;
    lambda_hat = lambda_t( tvec ,alpha_ref ,nhp_params.phi_n,nhp_params.nInd);
    fsize = 11;
    pcnt = 12;
    legarr={};
    %         pindArray = [4 , 8, 10 , numel(xsol)]; %batch2print
    pindArray = pcnt+[0, 1 , 2];%, 3]; %batch2print
    colArr = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.4660 0.6740 0.1880],	[0.4940 0.1840 0.5560]	};
    figure1 = figure(3);
    clf;
    tittext = sprintf('Visualizing updates conversion');
    title(tittext)
    hold all
    for p=1:numel(pindArray)
        legarr={};
        %             subplot(4,1,p), plot(tvec, Lambda_CVX, 'dk','MarkerSize',3,'MarkerIndices',1:4:length(tvec));
        subplot(numel(pindArray),1,p), plot(tvec, lambda_hat, '--k','LineWidth',1.5);
        hold on;
        ax=gca;
        ax.FontSize=fsize;
%         ax.TickDir = 'out'
        yticks([])
        legarr{1} = sprintf('$ \\lambda^*(t)$');
        pind = pindArray(p);
        subplot(numel(pindArray),1,p)
        hold on
        plot(tvec, lambda_t( tvec ,xbatchsol(:,pind+1),nhp_params.phi_n,nhp_params.nInd),...
            'Color',colArr{p},'LineWidth',2,'MarkerIndices',1:3:length(tvec),'LineStyle','-');
        %extract and plot events included in the current estimation
        It = [Tlim_div(1) Tlim_div(pind,2)];
        tK_inBatch = tK(tK>=It(1) & tK<=It(2));
        plot(tK_inBatch, 0.5*ones(size(tK_inBatch)), 'bd','MarkerFaceColor','y','MarkerSize',5,'Marker','d');
        xlabel('time[t]','FontSize',fsize)
%         ylabel('\lambda(t)')
        legarr{2} = sprintf('$ \\hat{\\lambda}_{%i}(t)$', (pind)*bksz);
        legarr{3} = 'Events';
        %print the legend
        h_leg =  legend(legarr,'Interpreter','latex');
        h_leg.FontSize=fsize;
        h_leg.Location='northeast';
        xlim([18 30])
        ylim([50 250])
%         xlim([0 10])
    end
end


% Create rectangle
annotation(figure1,'rectangle',...
    [0.30 0.72 0.18 0.16],...
    'LineStyle','-.');

% Create rectangle
annotation(figure1,'rectangle',...
    [0.30 0.42 0.18 0.16],...
    'LineStyle','-.');

% Create rectangle
annotation(figure1,'rectangle',...
    [0.30 0.12 0.18 0.16],...
    'LineStyle','-.');



% Add double arrows to where convergence occurs

% Create arrow
% annotation(figure1,'arrow',[0.498958333333333 0.498958333333333],...
%     [0.886662988966901 0.840521564694082]);
% 
% Create arrow
% annotation(figure1,'arrow',[0.498958333333333 0.498958333333333],...
%     [0.666001003009027 0.619859578736209]);
% 
% Create arrow
% annotation(figure1,'arrow',[0.5 0.5],[0.446342026078235 0.400200601805416]);
% 
% Create arrow
% annotation(figure1,'arrow',[0.498958333333333 0.498958333333333],...
%     [0.229692076228686 0.183550651955867]);
% 
% Create arrow
% annotation(figure1,'arrow',[0.498958333333333 0.498958333333333],...
%     [0.104312938816449 0.147442326980943]);
% 
% Create arrow
% annotation(figure1,'arrow',[0.498958333333333 0.498958333333333],...
%     [0.321965897693079 0.364092276830492]);
% 
% Create arrow
% annotation(figure1,'arrow',[0.498958333333333 0.498958333333333],...
%     [0.543630892678034 0.583751253761284]);
% 
% Create arrow
% annotation(figure1,'arrow',[0.5 0.499479166666667],...
%     [0.761283851554664 0.804413239719157]);
%% Plot Updates conversion - Column subfigures - Dynamic selection of time reference
% We go through the same plotting principal we did in the previous plot
% with the option to dynamically navigate through the time axis using the
% key inputs 'u/d' to move forward or backward.
if (1 || plotall)
    %%
    xbatchsol = XHIST_unlim;
    dirkey = 'u';
    pcnt = 1;
    lambda_hat = lambda_t( tvec ,alpha_ref ,nhp_params.phi_n,nhp_params.nInd);
    while(~strcmp(dirkey,'e'))
        legarr={};
        %         pindArray = [4 , 8, 10 , numel(xsol)]; %batch2print
        pindArray = pcnt+[0, 1 , 2, 3]; %batch2print
        colArr = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.4660 0.6740 0.1880],	[0.4940 0.1840 0.5560]	};
        figure(4),clf
        tittext = sprintf('Visualizing updates conversion');
        title(tittext)
        hold all
        for p=1:numel(pindArray)
            legarr={};
            %             subplot(4,1,p), plot(tvec, Lambda_CVX, 'dk','MarkerSize',3,'MarkerIndices',1:4:length(tvec));
            subplot(4,1,p), plot(tvec, lambda_hat, '-.k','LineWidth',1.5);
            %plotting the refrence estimate
            hold on;
            legarr{1} = sprintf('$ \\lambda^*(t)$');
            pind = pindArray(p);
            subplot(4,1,p),
            plot(tvec, lambda_t( tvec ,xbatchsol(:,pind+1) ,nhp_params.phi_n,nhp_params.nInd),...
                'Color',colArr{p},'LineWidth',1.5,'MarkerIndices',1:3:length(tvec));
            legarr{2} = sprintf('$ \\hat{\\lambda_{%i}(t)}$', pind*bksz);
            h_leg =  legend(legarr,'Interpreter','latex');
            h_leg.FontSize=12;
            h_leg.Location='northeast';
            %             xlim([0 30])
        end
        dirkey=input('Press f (Fwd) / b (Bwd) / e (Esc)\n','s');
        if (strcmp(dirkey, 'f'))
            pcnttmp = pcnt +1;
            pcnt = min(size(xbatchsol,2),pcnttmp);
        elseif(strcmp(dirkey, 'b'))
            pcnttmp = pcnt -1;
            pcnt = max(1, pcnttmp);
        end
    end
end

%% simulation for error vs buffer size
if (1 || plotall)
    bufErr_l = zeros(numel(ResArray{1}),K);
    for l=1:K
        bufErr_l(:,l) = ResArray{l} - xstrm_unlim;
    end
    %
    Y=1:110;
    X=1:10;
    Z=log10(bufErr_l.^2);
    figure(5),clf
    fh=figure(5)
    axes1 = axes('Parent',fh);
    [XX,YY]=meshgrid(X,Y);
    plot3(XX,YY,Z,'.-','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
    xh = xlabel('Max buffer size / block size')
    xh.FontSize = 15;
    yh = ylabel(sprintf('$ [t]$'),'Interpreter','latex');
    yh.FontSize=20;
    zh=zlabel(sprintf('$\\log \\| \\bar{x}_t-x^*_t \\|^2_2 $'),'Interpreter','latex')
    zh.FontSize=20;
    %     title('Truncation Error VS. Buffer Size')
    %     view(axes1,[12.5749999999999 23.9172413793104]); %3D view
    view(axes1,[0 0] ); % log vs buf size 2D view
    
    %%
    figure(6),clf
    fh=figure(6)
    axes6 = axes('Parent',fh);
    [XX,YY]=meshgrid(X,Y);
    plot3(XX,YY,Z,'.-','MarkerSize',15,'MarkerFaceColor','#D9FFFF')
    xh = xlabel('Max buffer size / block size')
    xh.FontSize = 15;
    yh = ylabel(sprintf('$ [t]$'),'Interpreter','latex');
    yh.FontSize=20;
    zh=zlabel(sprintf('$\\log \\| \\bar{x}_t-x^*_t \\|^2_2 $'),'Interpreter','latex')
    zh.FontSize=20;
    %     title('Truncation Error VS. Buffer Size')
    
    
    view(axes6,[90 0] ); % log vs buf size 2D view
    
    %%
    
    
end

if (0 || plotall)
    NMSE2 = sum(bufErr_l.^2,1)/size(bufErr_l,1);
    figure(7),clf
    semilogy(BUFSZARRAY, (NMSE2),'--bs','MarkerFaceColor','green');
    hx7 = xlabel('Buffer size');
    hx7.FontSize = 20;
    hy7 = ylabel(sprintf('$\\log {\\| \\underline{\\bar{x}}_T-\\underline{x}_T^* \\|^2_2} /{n(T+1)} $'),'Interpreter','latex');
    hy7.FontSize= 15;
    title('Truncation Error Vs Buffer Size')
end
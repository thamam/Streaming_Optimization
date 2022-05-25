%% Unified Results Ana;ysis
plotLineSections = @(T) plot(T.'*ones(1,100),linspace(0,1.1*params.Prate,100),'--k','linewidth',0.5);
colors=['r','g','b','m','k','--r','--g','--b','--m','--k',':r',':g',':b',':m',':k'];
SolList = solList_obj.list;
% SolList = {'FWD_W_single_BWD_UPDATEOut','CoordDescOut', 'WavesCoordDescOut','FWD__SWP__W__TAIL',...
%     'SWEEPOut','RANDOMCDOut', };
total = [];
solListLegend = strrep(SolList(2:end),'_','\_');
%% Assign solution output into analysis variables
for p=1:length(SolList)
    if strcmpi(SolList{p},'GlobalCvx')
        continue;
    end
    Out = eval(sprintf('Out_%s',SolList{p}));
    %         case 'CoordinateDescent'
    alphaValM = Out.alphaValM;
    W=size(Out.tK_sec_ind,1);
    K=Out.SectionIterationCount;
    activeFramesMat = Out.activeFramesMat;
    Out.Type = SolList{p};
    SolName = strrep( SolList{p},'_','\_');
    
    %% Print Activation Matrix
    if false
        %%
        figure,
        cMap=[0.2,0.2,.8;.9,.95,.2];
        imagesc(activeFramesMat)
        colormap(cMap)
        title(['Activation Matrix - ',SolName]);
        xlabel('Frames [w]');
        ylabel('Iterations')
    end
    
    if false
        %%  Print all results
        
        T=[params.nInd(1)-0.5+SIM{2}.options.sectionLength:SIM{2}.options.sectionLength:params.simTfinish];
        
        KK = 1:2:(K-1);
        %                 F= struct('cdata',[],'colormap',[]);
        %                 F(length(KK)) = struct('cdata',[],'colormap',[]);
        for k=KK
            activeFramesMask  = activeFramesMat(k,:);
            activeFramesList = find(activeFramesMask==1);
            figure(k), clf
            hold on
            plot(tt,Lambda_CVX.','--r','linewidth',1);
            alphaMat =  alphaValM(k+1,:);
            Lambda_k = lambda_t(tt,alphaMat.',params.phi_n,params.nInd);
             title(['Estimation Result - ',SolName]);
            h_legend = legend('\lambda_{Global}','\lambda_{k}')
            set(h_legend,'fontsize',24)
        end
    end
    
    if true
        %% Plot Final Results W/ global solution
        figure(3),clf
        plot(tt,Lambda_GroundTruth.','--g','linewidth',1.5);
        alphaMat =  alphaValM(K,:);
        Lambda_k = lambda_t(tt,alphaMat.',params.phi_n,params.nInd);
        hold all;
        plot(tt,Lambda_CVX.',':b','linewidth',1.5);
%         scatter(tK,Lambda_at_arrivalPoints,20,'*')
        scatter(tK,0*ones(size(tK)),20,'+r','LineWidth',1.5)
        Lambda_k = lambda_t(tt,alphaMat.',params.phi_n,params.nInd);
%         plot(tt,Lambda_k,'--r','linewidth',1);
        title(SolName) 
        %         plotLineSections(T)
        h_legend = legend('\lambda_{G. truth}','\lambda_{Global}','Arrival Times');%,'\lambda_{k}')
        set(h_legend,'fontsize',24)
    end
    
    
    %% Analysis number of times each frames is updated
    activeFramesMat = Out.activeFramesMat;
    frameActivationCount(p,:) = sum(activeFramesMat(1:K-1,:),1);
    %     h = figure(101),bar(1:W,h) ; hold all
    total(p) = sum(frameActivationCount(p,:));
    
end

%% Plot Histogram of frames activation
figure,
subplot(211),
framesInd = 1:W;
bar(framesInd,frameActivationCount(2:end,:).','grouped'); % groups by row
title('Frame Activation Histogram')
legend(solListLegend)

ax = subplot(212);
bar(total(2:end),'grouped'); % groups by row
ax.XTickLabel =solListLegend;

%% Analysis for JR Update

% comulative Error
alphaK = Out_FWD_STRM_BWD_UPDATE_JR.alphaValM(Out_FWD_STRM_BWD_UPDATE_JR.SectionIterationCount,:);
Lambda_JR = lambda_t(tt,alphaK.',params.phi_n,params.nInd);
cumErr = cumsum((Lambda_CVX(:).' -Lambda_JR(:).' ).^2./(1:numel(Lambda_JR)));
figure,
plot(tt,cumErr)
title('Comulative Error (normalized by num of samples')
xlabel('Time')
ylabel('e(t) = ( \int_0^t{( \lambda_{global-est}(t) - \lambda_{JR}(t))^2dt})^{1/2}')

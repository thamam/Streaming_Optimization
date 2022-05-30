%% Coordinate Descent plot of lambda function
%% Plotting .....
simOut = Out_wavesCD;
alphaValM = simOut.alphaValM; %CoordDescOut.alphaValM;
W=size(simOut.tK_sec_ind,1);
K=simOut.SectionIterationCount;
%%
close all
opts =  WCD_opts;
simOut.SectionIterationCount
T=params.simTstart:opts.sectionLength:params.simTfinish;
plotLineSections = @(T) plot(T.'*ones(1,100),linspace(0,1.1*params.Prate,100),'--k');
colors=['r','g','b','m','k','--r','--g','--b','--m','--k',':r',':g',':b',':m',':k'];
KK = [2,2+4:4:K-4,K];
if(~exist('saveVideo'))
    saveVideo = false
end

if saveVideo
    KK = 2:K;
else
    KK = [2,2+4:4:K-4,K];
end
%%
KK = 2:10:K;

for k=KK
    figure
    if exist('Lambda_CVX')
        plot(tt,Lambda_CVX.',':b','linewidth',1.5);
    end
    hold on;
    alphaMat =    alphaValM(:,k);
    Lambda_W    = lambda_t(tt,alphaMat.',params.phi_n,params.nInd);
    plot(tt,Lambda_W.','-r');
    plotLineSections(T)
    title(['Iteration k = ',num2str(k),'    Error Norm - ',num2str(CoordDescOut.alpha_diff_norm(k))]);
    legend('Global Sol','GCD SOl');
    F{k}=getframe(gcf); 
end
%
if saveVideo
    vidFileName = ['Vid_Res_',params.fileNhppName(1,1:end-3),'mp4'];
    v = VideoWriter(vidFileName);
    v.FrameRate = 1;
    open(v)
    for k=KK
        writeVideo(v,F{k})
    end
    close(v)
end


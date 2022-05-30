%%
[ttSup, T, plotLineSections, colors] = plot_basis_functions;

%%
W=4;
figure(7),clf
nIndSec=[1,5;6,10;11,15;16,20];

for w=1:W
    subplot(2,2,w)
%     ax=gca;    ax.XTick = params.nInd;
    hold all;
    for i=nIndSec(w,1):nIndSec(w,2)
        plot(ttSup,12*lambda_t(ttSup,1,params.phi_n,params.nInd(i)),colors(i))
    end
    plotLineSections(T)
    plot(ttSup,lambda_t(ttSup,alphaSyn,params.phi_n,params.nInd),'--k','linewidth',1)
    scatter(tK,Lambda_at_arrivalPoints,'*k');

%     scatter(tK,Lambda_at_arrivalPoints,'dk');

end
grid on
%%
plotPhi_nInd = @(i) plot(ttSup,lambda_t(ttSup,1,params.phi_n,params.nInd(i)),colors(i));

figure(8),clf
KW = [1:3,3:-1:1]
T=[2,6,10,14];

tSec=[1 5 ; 6 10; 11 15];
% for n=1:length(KW)
for n=1:3    
    w=KW(n);
%     subplot(numel(KW),1,n)
    subplot(3,1,n)
    ax=gca;    ax.XTick = params.nInd
    hold all;
    ylim([0 1])
    for i=tSec(w,1):tSec(w,2)
        plotPhi_nInd(i);
    end
    plotLineSections(T)
end
%%

function  plot_basis_functions(params)

L = params.splineOrder;
params.tsup;
nInd=params.nInd;
ttSup = linspace(params.tsup.Start,params.tsup.End,1000);
T=[params.nInd(1)-0.5:params.sectionLength:params.simTfinish];
plotLineSections = @(T) plot(T.'*ones(1,100),linspace(0,params.Prate,100),'--k');

colors=['r','g','b','m','k','--r','--g','--b','--m','--k',':r',':g',':b',':m',':k'];
figure(6),clf
ax=gca
hold on;
for i=1:params.nInd(end)
    plot(ttSup,params.Prate*lambda_t(ttSup,1,params.phi_n,params.nInd(i)),colors(1+mod(i,length(colors))))
end
plotLineSections(T)

grid 'on'
ax.XTick = params.nInd;
end
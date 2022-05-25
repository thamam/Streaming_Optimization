function [Fi,gmap,Zgind,InitValues,xupdateGrad] =  prepadmmNHPP(P,params,tK)
%% documantation


%% Initialize
% rho = params.rho;
L=params.splineOrder;
Ti = settimeintervals(params.simTstart ,params.simTfinish,P);
%% Divide samples/Set global variable map/Set local objective functions
for i=1:P
    % Generate gloval variable mapping and subsets w.r.t. time intervals
    % assigned to each agent
    tempZg = [Ti{i}(1)-(L+1)/2:1:Ti{i}(2)+(L+1)/2];
    Zgind{i} = tempZg(2:end-1); %psi indexes corresponding to global variable order.i.e. (x_i)_j=zgind(zgmap(i,j))
    [~,gmap{i}]=ismember(Zgind{i},params.nInd) ;%global var map
    bi{i} = tK(tK>=Ti{i}(1) & tK<=Ti{i}(2)); %divide samples
    InitValues.x0{i} = ones(numel(gmap{i}),1) ;
    InitValues.y0{i} = ones(numel(gmap{i}),1) ;
    InitValues.z0{i} = ones(numel(gmap{i}),1) ;   
    ZgL{i} = length(gmap{i});
    [ ci,PSIi ] = basisFunctionPrep(Ti{i},bi{i},Zgind{i},params.phi_n); %generate psi dependent constants
    Fi{i} = @(x,yik,zetaik,rho) localneglogObjectiveFunc(x,yik,zetaik,rho,ci,PSIi)  ;
    xupdateGrad{i} = @(x,yik,zik,rho) ci-    sum(PSIi*diag(1./(x.'*PSIi)),2) + yik + rho*(x-zik);
end



    





end

%%%% -----------------------------------------------------------------------------------------------------------------------------%%%%
%%%% -----------------------------------------------------------------------------------------------------------------------------%%%%
%%%%                                                                                 Subroutines                                                                              %%%%
%%%% -----------------------------------------------------------------------------------------------------------------------------%%%%
%%%% -----------------------------------------------------------------------------------------------------------------------------%%%%
%%
function [ ci,PSIi ] = basisFunctionPrep(Ti,tK_i,Zgi,phi_n)
%% Doc:
% c{i} := integral_{t\in T_i} {Psi_i}_j(t)dt   
    nInd=Zgi;%
    d= numel(nInd);
    K = numel(tK_i);
    ci = zeros(d,1);
    for i=1:d %compute integral consts
        ci(i) = integral(@(t)phi_n(t,nInd(i)),Ti(1),Ti(2));
    end
    
%Assemble PSI | PSI[i,j] = psi_n(tK[j])
    PSIi = zeros( d,K);
    for i=1:numel(nInd)
        PSIi(i,:) = phi_n(tK_i,nInd(i));
    end
end

function [val] = localneglogObjectiveFunc(x,yik,zetaik,rho,ci,PSIi)
    val = ( x.'*ci - sum(log(x.' * PSIi)) +yik.' * x + (rho/2)*norm(x-zetaik,2)^2); %local objective function
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ti] = settimeintervals(Tstart,Tfinish,P)
    sampleIntervalsticks = round(linspace(Tstart,Tfinish,P+1)); %samples Intervals
    tempTicks = sort(repmat(sampleIntervalsticks,[1,2]));
    tempTicks = tempTicks(2:end-1);
    TiMat=reshape(tempTicks,2,[]).';
    Ti = mat2cell(TiMat,ones(1,P),2);
end
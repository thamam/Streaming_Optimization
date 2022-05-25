function [xsol, it_hist, outstat,dbgdata] =  nhpp_admmconsensus(F,InitValues,parms,gmap,Zgind,h_Xupd,h_Zupd,xupdateGrad)
% Documentation
% f_i= F{i}; cell of pointers to the local objective functions
% x0 = cell of size P x 1, where the i-th element of the cell x^0_i , i.e. the initial value for the x_i

%% Initialize Code

% initialize it_hist, ierr, and set the iteration parameters
% maxit = 20; rtol=1e-4; ierr = 0;
N = parms(4);
maxit = parms(1) ; rtol = parms(2) ;
rho = parms(3);
it_hist=[];
k=0; %
P = numel(F);
cnsnorm = @(z,x,y)  norm( z-h_Zupd(x,y,N,P,gmap,rho ),2);
atol=rtol;
x0 = InitValues.x0;
zeta = InitValues.z0;
y=InitValues.y0;
x=x0;
z=zeros(N,1);

fnrm0= cnsnorm(z(:,1),x(1,:),y(1,:)) ;
fnrm=fnrm0;
it_hist(k+1)=fnrm;
outstat(k+1, :) = [k fnrm fnrm/fnrm0];
stop_tol=atol;% + rtol*fnrm;

%% Init

fmin_options = optimoptions(@fminunc,'Algorithm','trust-region',...
    'SpecifyObjectiveGradient',true);
while(fnrm/fnrm0 > stop_tol && k < maxit)
    % The stopping criteria used is either maxitr or norm of the consensus
    % residual among th agnets. Obviously, if they reach full consensus the
    % algorithm can terminate.
    k=k+1;
    %% In parallel - Compute Xi
    for i=1:P % update prime variable x
        %Fi{i} @(rho,zetaik,yik)   -Fi header
        fi = @(x) F{i}(x,y{k,i},zeta{k,i},rho);  %function val handle
        gi = @(x) xupdateGrad{i}(x,y{k,i},zeta{k,i},rho);  %gradient handle
        x{k+1,i}  = fminunc(@(x) localobjectiveinfo(x,fi,gi) ,x{k,i},fmin_options);
    end
    
    %% Communicate
    % set the weight for the averaging - doesn't require any prior
    % knowledge how data is partitioned
    
    z(:,k+1) = ProxProj( h_Zupd(x(k+1,:),y(k,:),N,P,gmap,rho ));
    
    % map Z to zeta
    for i=1:P
        zeta{k+1,i} = z(gmap{1,i},k+1);
    end
    
    %% Aggregate
    for i=1:P % Update dual variable u
        y{k+1,i} = y{k,i} + rho*(x{k+1,i}-zeta{k+1,i} );
    end
    
    %     if k==1 %set the zero norm as the norm of the first result.
    %         fnrm0=norm(z(:,k+1));
    %     end
    fnrm = cnsnorm(z(:,k),x(k+1,:),y(k,:)) ; %notice z here is array and not cell  !!!
    it_hist(k+1)=fnrm;
    outstat(k+1, :) = [k fnrm fnrm/fnrm0];
end
% Comprise final solution
xsol = z(:,k+1);  % z contains the closest version to the consensus the alg. has reached.
%% Assign debug data
dbgdata.z = z;
dbgdata.x=x;
dbgdata.y=y;
dbgdata.gmap = gmap;

end



%%%% -----------------------------------------------------------------------------------------------------------------------------%%%%
%%%% -----------------------------------------------------------------------------------------------------------------------------%%%%
%%%%                                                                                 Subroutines                                                                              %%%%
%%%% -----------------------------------------------------------------------------------------------------------------------------%%%%
%%%% -----------------------------------------------------------------------------------------------------------------------------%%%%
%%
function [f,g] = localobjectiveinfo(x,fi,gi)
% Calculate objective f
f = fi(x);
if nargout > 1 % gradient required
    g = gi(x);
end
end

function [ zplus] = ProxProj(z)
zplus = z;
zplus(zplus<0)=0;
end
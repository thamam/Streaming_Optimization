function [varargout] = NWSolver(Jobj,x, varargin)
% NWSolver.m 
% [sol] =  NWSolver(Jobj,x,Ts,nonegxMod,t, mu, eps)
% Return the Logbarrier constrained Newton solver

% mest is an estimate of the number of constraints-good enough to compute the stopping criteria
% mest = round(diff(Jobj.support)/Ts);

%log barrier default values
t0 = 1; mu = 10; eps=1e-4;

dbgstr=struct;
bksz = Jobj.bksz;
NT = Jobj.NT;
T=Jobj.T;
MAXBUFSIZE = Jobj.bfsz;
Tmbp1 = max(1,MAXBUFSIZE+2-T);
supp = Jobj.support;


%Initialize

% Determine if trucation is active
trncEn = (T>MAXBUFSIZE);

if trncEn  %truncation is active, f_tild is used
    xf_Tmb = x(1:bksz);
    f= @(y) objJ_fungrad(y,Jobj,trncEn, xf_Tmb);
    x0=x(bksz+1:end);
else
    f= @(x) objJ_fungrad(x,Jobj, trncEn);
    x0=x;
end


%Active basis func indices
xInd = (Jobj.xIndtbl(1):1:Jobj.xIndtbl(end)).';


%% Prepare Newton Solver
% addpath(genpath('../../../CODES/'));
% % Input
isham = 1 ;  rsham = 0 ; % is Newton's method
maxit = 80; % maxmium number of nonlinear iterations
atol = 1.d-6; rtol = 1.d-6;
tol = [atol, rtol] ; %relative/absolute tolerance

parms = [maxit, isham, rsham];

options.UseHessian=1;
options.bksz = bksz;

%% Log barrier initialization
if nargin >3
    t0 = varargin{1};
    if nargin >4
        mu = varargin{2};
        if nargin > 5
            eps = varargin{3};
        end
    end
end

%% Solver - If it's the first iteration
if T==2
    % p^* is the value of the lagrangian at x^*
    % f0 is value of original objective at x^* (can be approximated by
    % evaluating f at very large t)
    % Stopping criterion:
    %    (1) if m/t < eps ; m=dim(phi), or
    %    (2) if (f0(x) - pstar)<= eps
    
    
    [y,iterHist] = NwBarrier(x, f, parms, options, tol,t0, mu, eps);
    
    
else  % T >=3
    
    % Finding y0_T
    fT =  Jobj.fObjArray.datalist{end};
    y0_Tm1 = x0(end-2*bksz+1:end-bksz);
    y0_T_tmp = x0(end-bksz+1:end);
    %     w0 = x0(end-2*bksz+1:end);
    %     Aeq_T = [eye(bksz),zeros(bksz)]; beq_T = y0_Tm1;
    %
    %     fTtild = @(y) objft_fungrad(y,fT);
%     fTtild_ = @(yT) objft_parder([y0_Tm1;yT],fT);
    
%     [y0_T] = NwBarrier(y0_T_tmp, fTtild_, parms, options, tol,t0, mu, eps);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    A_T  = -eye(bksz*2);

    b_T = zeros(size(A_T,1),1);
   
    Aeq_T = [eye(bksz),zeros(bksz)]; beq_T = y0_Tm1;
    
    w0 = x0(end-2*bksz+1:end);
    
    fT_h = @(y) objft_fun(y,fT);
    
    yT_mp = fmincon(fT_h, w0,A_T,b_T,Aeq_T,beq_T,[],[],[]);
    
    y0_T = yT_mp(bksz+1:end);
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Solve for y*
    % form y0
    y0=x0; y0(end-bksz+1:end) = y0_T;
    %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    
    [y, iterHist] = NwBarrier(y0, f, parms, options, tol,t0, mu, eps);      
    
    %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    
end
ystar = y;

varargout{1}=ystar;

if nargout>1
    varargout{2} = iterHist;
    if nargout>2
    [Fystar, Hystar] = f(ystar);
    varargout{3} = Fystar;
    if nargout>3
        varargout{4} = Hystar ;
    end
    end

end






end

% function [gradf, hesf] = objJ_fungrad(x, fobj, trunEn, xf_Tmb)
% 
% if trunEn %we modify the function call by fixing the value of xf_Tmb
%     N = numel(x);
%     bksz = fobj.bksz;
%     x_act_ind = bksz+(1:N);
%     
%     xtld = [xf_Tmb; x];
% 
%     tmp_grad = fobj.GradJ(xtld);
%     gradf = tmp_grad(x_act_ind);
% 
%      if nargout  > 1
%         tmp_hesf = fobj.HessJ(xtld);
%         hesf = tmp_hesf(x_act_ind,x_act_ind);
%      end
% else
%     gradf = fobj.GradJ(x);
%     % Gradient of the objective function:
%     if nargout  > 1
%         hesf = fobj.HessJ(x);
%     end
% end
% 
% end

% function [gradf, hesf] = objft_fungrad(x, fobj)
% gradf = fobj.fgrad(x);
% % gradf = gradf_(fobj.ntm1+1:end);
% % Gradient of the objective function:
% if nargout  > 1
%     hesf = fobj.fhessian(x);
%     %     hesf = hesf_(fobj.ntm1+1:end,fobj.ntm1+1:end);
% end
% 
% end

% function [gradf_xt,hes_xt] = objft_parder(x,fobj)
% 
% N = fobj.ntm1+fobj.nt;
% 
% U = eye(N);
% Ut = U(fobj.ntm1+1:end,:);
% 
% gradf = fobj.fgrad(x);
% gradf_xt = Ut*gradf;
% 
% if nargout  > 1
%     hesf = fobj.fhessian(x);
%     hes_xt = Ut*hesf*Ut.';
% end
% 
% 
% end

% function [Jobj]=updatearrayprop(Jobj, pval, Tmbp1, MAXBUFSIZE)
% %update property value over a cellarray of class
% % obj is handle to the class object cell array
% % pname is the name of the property field to be upodated
% % pval is the new value to be assigned
% for i=Tmbp1:MAXBUFSIZE
%     propmodcom = sprintf('Jobj.fObjArray.datalist{%i}.ltObj.lb_delta=pval ;',i);
%     eval(propmodcom);
% end
% end

% function [f,gradf] = objft_fun(x, fobj)
% f = fobj.feval(x);
% % gradf = gradf_(fobj.ntm1+1:end);
% % Gradient of the objective function:
% if nargout  > 1
%     gradf = fobj.fgrad(x);    %     hesf = hesf_(fobj.ntm1+1:end,fobj.ntm1+1:end);
% end
% 
% end

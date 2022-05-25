function [varargout]= NwBatchSolver(x0, f_batch, xInd,...
    supp, Ts,options,varargin)
%Solve the conditiondional batch NHPP optimization program from 0 to any
%give T\leq Tend

mest = round(diff(supp)/Ts);

%log barrier default values
t0 = 1; mu = 10; eps=1e-10;

if nargin>6
    t0 = varargin{1};
    if nargin >7
        mu=  varargin{2};
        if nargin>8
            eps= varargin{3};
        end
    end
end


% nonegxMod = true;
dbgstr =struct;

Tgrid=[];
% Prep objective
Obj = f_batch.ltObj;
Tlim = Obj.Tlim;
tK = Obj.tK_t;
% minimize x∈X〈x,a〉−∑_m log (〈x,cm〉)
NT= numel(xInd);

%%LGBR settings
Tgrid = (supp(1)+Ts):Ts:(supp(2)-Ts);
At = -f_batch.ltObj.step(Tgrid,xInd(:));


a = zeros(NT,1);
for j=1:NT
    a(j) = integral(@(t) Obj.step(t,xInd(j)),Tlim(1),Tlim(2));
end

cm = Obj.step(tK,xInd);

% feval = @(x) a.'*x  -sum ( log(cm.'*x)) ;
fgradobj = @(x) a -  sum(cm ./ (cm.'*x).',2);
lb_grad_f = @(x)  -sum(At ./ ((At.'*x).'),2);
fgrad = @(x,t) t*fgradobj(x) + lb_grad_f(x);

fhess = @(x,t) HessianComp(x,cm,At,t);

%% Prepare Newton Solver
% addpath(genpath('../../../CODES/'));
% % Input
isham = 1 ;  rsham = 0 ; % is Newton's method
maxit = 80; % maxmium number of nonlinear iterations
atol = 1.d-6; rtol = 1.d-6;
tol = [atol, rtol] ; %relative/absolute tolerance

parms = [maxit, isham, rsham];

options.UseHessian=1;


% Call for fmincon solver
t = t0/mu; %init t W/ precompensate for the first iteration.
lgbrCnt = 0;
y0=x0;

while ( (mest/t)>eps || lgbrCnt<1 ) %making sure we solve for x at least one time
    % or m/t < eps
    
    % 0. Update counter
    lgbrCnt  = lgbrCnt +1;
    
    % 1. Set/Increase t:  t := µt
    t = mu * t;
   
    % 2. compute x* =argmin t*f + phi
    f = @(x) objJ_fungrad(x, t, fgrad, fhess);

    [y, it_hist, ierr,outstat] = nsold(y0, f,tol,parms, options);
    
    % 3. Update. x := x*(t).
    y0 = y;
end



ystar = y;

varargout{1}=ystar;

if nargout>1
    [Fystar, Hystar] = f(ystar);
    varargout{2} = Fystar;
    if nargout>2
        varargout{3} = Hystar ;
    end
end


end



%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%
%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%
%%%                         Auxilary functions
%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%
%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%
function [gradf, hesf] = objJ_fungrad(x,t, gradf_h, hes_h)

gradf = gradf_h(x,t);
% Gradient of the objective function:
if nargout  > 1
    hesf = hes_h(x,t);
end
end

function [HesOut] =HessianComp(x,cm,At,t)


d = numel(x);
if nargin==1
    x = obj.x;
end

%computing the objective Hessian
AtT=cm.';
den2 = ((AtT*x).').^2 ; % den = [ sum(x.'*p0) , sum(x.'*p_1),...,sum(x.'*p_K)]
Hout=zeros(d);
for i=1:d
    a_ti = AtT(:,i).' ;
    for j=1:d
        a_tj = AtT(:,j).' ;
        Hout(i,j) = sum(a_ti.*a_tj./den2);
    end
end


%computing the LGBR Hessian
%computing the objective Hessian
AtT=At.';
den2 = ((AtT*x).').^2 ; % den = [ sum(x.'*p0) , sum(x.'*p_1),...,sum(x.'*p_K)]
Hlgb=zeros(d);
for i=1:d
    a_ti = AtT(:,i).' ;
    for j=1:d
        a_tj = AtT(:,j).' ;
        Hlgb(i,j) = sum(a_ti.*a_tj./den2);
    end
end

%% combine the Hessians
HesOut = t*Hout + Hlgb;


end
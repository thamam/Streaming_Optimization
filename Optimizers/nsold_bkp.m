function [sol, it_hist, ierr,outstat] = nsold_bkp(x,f,tol,parms, options)
% This is modification of nsol.m to support analytic Hessian
% Newton solver, locally convergent
% solver for f(x) = 0
%
% Hybrid of Newton, Shamanskii, Chord
%
% C. T. Kelley, November 26, 1993
%
% This code comes with no guarantee or warranty of any kind.
%
% function [sol, it_hist, ierr] = nsol(x,f,tol,parms)
%
% inputs:
%        initial iterate = x
%		 function = f
%        tol = [atol, rtol] relative/absolute
%			error tolerances
%	 parms = [maxit, isham, rsham]
%			maxit = maxmium number of iterations
%				default = 40
%		isham, rsham: The Jacobian matrix is
% 		computed and factored after isham
%               updates of x or whenever the ratio
%		of successive infinity norms of the
%               nonlinear residual exceeds rsham.
%			isham = 1, rsham = 0 is Newton's method,
%			isham = -1, rsham = 1 is the chord method,
% 			isham = m, rsham = 1 is the Shamanskii method
%      	       defaults = [40, 1000, .5]
%
% output:
%	sol = solution
%	it_hist = infinity norms of nonlinear residuals
%			for the iteration
%	ierr = 0 upon successful termination
%	ierr = 1 if either after maxit iterations
%             the termination criterion is not satsified
%             or the ratio of successive nonlinear residuals
%             exceeds 1. In this latter case, the iteration
%				is terminted.
%
%
% internal parameter:
%       debug = turns on/off iteration statistics display as
%               the iteration progresses
%
% Requires: diffjac.m, dirder.m
%
% Here is an example. The example computes pi as a root of sin(x)
% with Newton's method and plots the iteration history.
%
%
%  x=3; tol=[1.d-6, 1.d-6]; params=[40, 1, 0];
%  [result, errs, it_hist] = nsol(x, 'sin', tol, params);
%  result
%  semilogy(errs)
%

%
% set the debug parameter, 1 turns display on, otherwise off
%
debug=0;
% 
% Global  line search parameters set alpha, sigma0, sigma1, maxarm, and restart_limit
%
alpha = 1.d-4; sigma0=.1; sigma1=.5; maxarm = 50; 

%
% initialize it_hist, ierr, and set the iteration parameters
%
ierr = 0;
maxit=40;
isham=1000;
rsham=.5;
if nargin == 4
    maxit=parms(1); isham=parms(2); rsham=parms(3);
end
rtol=tol(2); atol=tol(1);
% it_hist=[];
itc=0;

n = length(x);
fnrm=1;

it_histx(itc+1,1)=fnrm; it_histx(itc+1,2)=0; it_histx(itc+1,3)=0;
outstat(itc+1, :) = [itc fnrm 0 0 ];


%
% evaluate f at the initial iterate
% compute the stop tolerance
%
if options.UseHessian
    [f0, H0] = feval(f,x);
else
    [f0] = feval(f,x);
end
% it_hist=[it_hist,fnrm];
fnrmo=1;
fnrm=fnrmo ;

itsham=isham;
stop_tol=atol+rtol*fnrm;
%
% main iteration loop
%
while(fnrm > stop_tol & itc < maxit)
    %
    % keep track of the ratio (rat = fnrm/frnmo)
    % of successive residual norms and
    % the iteration counter (itc)
    %
    rat=fnrm/fnrmo;
%     outstat(itc+1, :) = [itc fnrm rat];
    fnrmo=fnrm;
    itc=itc+1;
    %
    % evaluate and factor the Jacobian
    % on the first iteration, every isham iterates, or
    % if the ratio of successive residual norm is too large
    %
    if(itc == 1 | rat > rsham | itsham == 0)
        itsham=isham;
        if options.UseHessian
            [l,u] =   lu(H0);
        else
            [l, u] = diffjac(x,f,f0);
        end
    end
    itsham=itsham-1;
    %
    % compute the step
    %
    tmp = -l\f0;
    step = u\tmp;
    %% Global convergence - line search
    %
    %   The line search starts here.
    %
    xold=x;
    lambda=1; lamm=1; lamc=lambda; iarm=0;
    xt = x + lambda*step;
    ft=feval(f,xt);
    nft=norm(ft); nf0=norm(f0); ff0=nf0*nf0; ffc=nft*nft; ffm=nft*nft;
    while nft >= (1 - alpha*lambda) * nf0
        %
        %   apply the three point parabolic model
        %
        if iarm == 0
            lambda=sigma1*lambda;
        else
            lambda=parab3p(lamc, lamm, ff0, ffc, ffm);
        end
        %
        % update x; keep the books on lambda
        %
        xt=x+lambda*step;
        lamm=lamc;
        lamc=lambda;
        %
        % keep the books on the function norms
        %
        ft=feval(f,xt);
        nft=norm(ft);
        ffm=ffc;
        ffc=nft*nft;
        iarm=iarm+1;
        if iarm > maxarm
            disp(' Armijo failure, too many reductions ');
            ierr=2;
            disp(outstat)
            it_hist=it_histx(1:itc+1,:);
            sol=xold;
            return;
        end
    end
    x=xt;
    %Reevaluate to get both Gradient and Hessian
    if options.UseHessian
        [f0, H0] = feval(f,x);
    else
        [f0] = feval(f,x);
    end    %
    %   end of line search
    %          
    
    fnrm=norm(f0,inf);
    it_histx(itc+1,1)=fnrm; 
    
%    
%   How many function evaluations did this iteration require?
%
    it_histx(itc+1,2)=it_histx(itc,2)+iarm+1;
    if(itc == 1) it_histx(itc+1,2) = it_histx(itc+1,2)+1; end
    it_histx(itc+1,3)=iarm;
    
    
    rat=fnrm/fnrmo;
    if debug==1
        disp([itc fnrm rat])
        %         debugRes{itc}.rat = rat;
        %         debugRes{itc}. = rat;
    end
    outstat(itc+1, :) = [itc fnrm rat iarm];
    %
    % if residual norms increase, terminate, set error flag
    % %
    %     if rat >= 1
    %         ierr=1;
    %         sol=xold;
    %         disp('increase in residual')
    %         disp(outstat)
    %         return;
    %     end
    % end while
end
sol=x;
it_hist=it_histx(1:itc+1,:);

if debug==1
    disp(outstat)
    it_hist=it_histx(1:itc+1,:);
end
%
% on failure, set the error flag
%
if fnrm > stop_tol
    ierr = 1;
end

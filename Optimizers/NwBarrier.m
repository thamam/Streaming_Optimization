function [yout, it_hist] = NwBarrier(y0, f0, parms, options, tol,varargin)
%NwBarrier Summary of this function goes here
%   Detailed explanation goes here
% varargin = {t0, mu, eps}


t0 = 1; mu = 10; eps=1e-10;

m = numel(y0);

if nargin >5
    t0 = varargin{1};
    if nargin >6
        mu = varargin{2};
        if nargin > 7
            eps = varargin{3};
        end
    end
end

t = t0/mu; %init t W/ precompensate for the first iteration.
lgbrCnt = 0;
frat = m/t;

while ( frat>eps || lgbrCnt<1 ) %making sure we solve for x at least one time
    % or m/t < eps
    
    % 0. Update counter
    lgbrCnt  = lgbrCnt +1;
    
    % 1. Set/Increase t:  t := µt
    t = mu * t;
    
    % 2. Update objective with barrier weight
    f = @(x) lgbrObj_Comp(x,f0, t);
    
    % 3. compute x* =argmin t*f + phi
    [sol_y,it_hist, ierr,outstat ] = nsold(y0, f ,tol ,parms, options);
    
    % 4. Update. x := x*(t).
    y0 = sol_y;
    
    frat = m/t;
end

yout  = y0;
end



function [gradf, hesf] = lgbrObj_Comp(x,f0, t)
%[f, gradf, hes_f] = lgbrObj_Comp(x,f0_obj, N) returns the value,gradient,
%and Hessian of f_0 with the addition of the log barrier function
%Input -
%  t is the log barrier weight
% f0 returns [grad, Hes] based on number of outputs
%Output
% f(x) = t*f0(x) + phi(x) where phi(x)=-sum_m log(f_m) is the log barrier
% function and f_m is the constraint function We use here f_m<=0 => f_m(x)
% = -x_m , wich we use to impose the constraint x_m>=0

% phi_x = @(x) -sum(log(x));
gradphi = @(x) -1./x(:);

if nargout == 1
    gradf = t*f0(x)+ gradphi(x); 
else
    [gradf0, hesf0]=f0(x);
    gradf = t*gradf0 + gradphi(x); 
    
    hesphi = @(x) diag(1./(x.^2));
    hesf = t* hesf0 + hesphi(x) ;
end

end


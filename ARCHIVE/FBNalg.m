%  Level 2
function [xstar,outputArg2] = FBNalg(y0,fObj,parms)
%RTNsolver Summary of this function goes here
%   Detailed explanation goes here

nT = fObj.nT ; % total # of variables

%% Initialize 
if numel(y0)<nT % check if y_T\nw{0} is given, if not compute it 
   y0_T = argmin(@w f_T(y0_Tm1, w));
   y=[y0(:);y0_T(:)];
else
    y=y0;
end

%% Newton solver for remaining of code

[xstar, it_hist, ierr] = nsola(y,F,tol,parms);



end
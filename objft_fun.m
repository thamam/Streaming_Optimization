function [f,gradf] = objft_fun(x, fobj)
f = fobj.feval(x);
% gradf = gradf_(fobj.ntm1+1:end);
% Gradient of the objective function:
if nargout  > 1
    gradf = fobj.fgrad(x);    %     hesf = hesf_(fobj.ntm1+1:end,fobj.ntm1+1:end);
end

end
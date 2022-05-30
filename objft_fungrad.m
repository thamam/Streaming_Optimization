function [gradf, hesf] = objft_fungrad(x, fobj)
gradf = fobj.fgrad(x);
% gradf = gradf_(fobj.ntm1+1:end);
% Gradient of the objective function:
if nargout  > 1
    hesf = fobj.fhessian(x);
    %     hesf = hesf_(fobj.ntm1+1:end,fobj.ntm1+1:end);
end

end

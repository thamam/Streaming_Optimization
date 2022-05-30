function [gradf_xt,hes_xt] = objft_parder(x,fobj)

N = fobj.ntm1+fobj.nt;

U = eye(N);
Ut = U(fobj.ntm1+1:end,:);

gradf = fobj.fgrad(x);
gradf_xt = Ut*gradf;

if nargout  > 1
    hesf = fobj.fhessian(x);
    hes_xt = Ut*hesf*Ut.';
end


end
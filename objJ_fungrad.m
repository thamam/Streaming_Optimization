function [gradf, hesf] = objJ_fungrad(x, fobj, trunEn, xf_Tmb)

if trunEn %we modify the function call by fixing the value of xf_Tmb
    N = numel(x);
    bksz = fobj.bksz;
    x_act_ind = bksz+(1:N);
    
    xtld = [xf_Tmb; x];

    tmp_grad = fobj.fgrad(xtld);
    gradf = tmp_grad(x_act_ind);

     if nargout  > 1
        tmp_hesf = fobj.fhessian(xtld);
        hesf = tmp_hesf(x_act_ind,x_act_ind);
     end
else
    gradf = fobj.fgrad(x);
    % Gradient of the objective function:
    if nargout  > 1
        hesf = fobj.fhessian(x);
    end
end

end
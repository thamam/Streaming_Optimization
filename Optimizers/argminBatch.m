function [xstar, fval,dbgstr] = argminBatch(x0, ...
    f_batch, xInd, supp, Ts,options, nonegxMod)
%Solve the conditiondional batch NHPP optimization program from 0 to any
%give T\leq Tend

% nonegxMod = true;
dbgstr =struct;

Tgrid=[];
% Prep objective
Obj = f_batch.ltObj;
Tlim = Obj.Tlim;
tK = Obj.tK_t;
% minimize x∈X〈x,a〉−∑_m log (〈x,cm〉)
NT= numel(xInd);

a = zeros(NT,1);
for j=1:NT
    a(j) = integral(@(t) Obj.step(t,xInd(j)),Tlim(1),Tlim(2));
end

cm = Obj.step(tK,xInd);

feval = @(x) a.'*x  -sum ( log(cm.'*x)) ;
fgrad = @(x) a -  sum(cm ./ (cm.'*x).',2);
fhess = @(x) HessianComp(cm, x);
fhandle = @(x) objJ_fungrad(x, feval, fgrad, fhess);


%% prepare inequlaity constaints
if(nonegxMod)
    A=-eye(numel(x0));
    b = zeros(size(A,1),1);
else
    %prepare inequality constraints
    Tgrid = (supp(1)+Ts):Ts:(supp(2)-Ts);
    A = -f_batch.ltObj.step(Tgrid,xInd).';
    b = zeros(size(A,1),1);
end

% Call for fmincon solver
[xstar, fval, exitflag, output, lambda, grad, hessian] = ...
    fmincon(fhandle ,x0,A,b,[],[],[],[],[],options);


dbgMod=true;

if (dbgMod) %
    dbgstr.A=A;
    dbgstr.x0=x0;
    dbgstr.a=a;
    dbgstr.cm=cm;
    dbgstr.Tgrid=Tgrid;
end



end

function [f,gradf, hesf] = objJ_fungrad(x, f_h,gradf_h, hes_h)
f = f_h(x);
% Gradient of the objective function:
if nargout  > 1
    gradf = gradf_h(x);
    if nargout > 2
        hesf = hes_h(x);
    end
    
end
end

function [HesOut] =HessianComp(cm,x)


d = numel(x);
if nargin==1
    x = obj.x;
end
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

HesOut = Hout;

end
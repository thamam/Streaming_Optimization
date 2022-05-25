function [xstar, fval,y0, dbgstr ] = argmincon(Jobj,x,Ts,fminConoptions, options)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% If nonegxMod==true, use inequality constraints to force x s.t. x>=0. If
% false, then we use the
%ineq constraints that x.'^psi(t)>=0 over a finite grid.
% nonegxMod = true;


%% Initialize
A=[]; b=[]; Aeq=[]; beq=[];Tgrid=[];

nonegxMod = options.nonegxMod;
dbgstr=struct;
bksz = Jobj.bksz;
NT = Jobj.NT; 
T=Jobj.T;
MAXBUFSIZE = Jobj.bfsz;
Tmbp1 = max(1,MAXBUFSIZE+2-T);
supp = Jobj.support;




% Determine if trucation is active
trncEn = (T>MAXBUFSIZE);

if fminConoptions.SpecifyObjectiveGradient
    f = @(x) objJ_fungrad(x,Jobj); 
else
    f = @(x) Jobj.J_Eval(x);
end
x0=x;

%Actiuve basis func indices
xInd = (Jobj.xIndtbl(1):1:Jobj.xIndtbl(end)).';

if T==2
    y0=x0; %for the first round only - no pre-solver
    f2 = Jobj.fObjArray.datalist{end};
    
    %prepare inequality constraints
    if(nonegxMod)
        A=-eye(numel(y0));
        b = zeros(size(A,1),1);
    else
        Tgrid = (supp(1)+Ts):Ts:(supp(2)-Ts);
        A = -f2.ltObj.step(Tgrid,xInd).';
    end
    b = zeros(size(A,1),1);
    
    [xstar,fval] = fmincon(f,y0,A,b,[],[],[],[],[],fminConoptions);

    
else % T >=3
    % Finding y0_T
    fT =    Jobj.fObjArray.datalist{end};
   
    if(nonegxMod)
        A_T  = -eye(bksz*2);
    else
        xInd_T = [fT.ltObj.xtm1_ind,fT.ltObj.xt_ind];
        supp_T = fT.ltObj.f_sup;
        Tgrid_T = (supp_T(1)+Ts):Ts:(supp_T(2)-Ts);
        A_T = -fT.ltObj.step(Tgrid_T,xInd_T).';
    end
    b_T = zeros(size(A_T,1),1);
    
    y0_Tm1 = x0(end-2*bksz+1:end-bksz);
    
    Aeq_T = [eye(bksz),zeros(bksz)]; beq_T = y0_Tm1;
    
    w0 = x0(end-2*bksz+1:end);
    
    fT_h = @(y) objft_fungrad(y,fT);
    
    y_tmp = fmincon(fT_h, w0,A_T,b_T,Aeq_T,beq_T,[],[],[],fminConoptions);
    
    y0=x0; y0(end-bksz+1:end) = y_tmp(end-bksz+1:end);
    % --------------------------------------------------------------------
    % Now we are ready to solve for y*
    % --------------------------------------------------------------------
    %% prepare equality constraints
    if trncEn %(only when trucation is active)
        Aeq = [eye(bksz), zeros(bksz, bksz*(min(T,MAXBUFSIZE+1)-1))];
        beq = y0(1:bksz);
    end
    
    %prepare inequality constraints
    if(nonegxMod)
        A=-eye(numel(y0));
        b = zeros(size(A,1),1);
    else
        Tgrid = (supp(1)+Ts):Ts:(supp(2)-Ts);
        A = -fT.ltObj.step(Tgrid,xInd).';
        b = zeros(size(A,1),1);
    end
    
    %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    [xstar,fval,exitflag,output,lambda,grad,hessian] = ...
        fmincon(f,y0,A,b,Aeq,beq,[],[],[],fminConoptions);
    %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
end

dbgmod= false;
if (dbgmod)
    dbgstr.y0=y0;
    dbgstr.A=A;
    dbgstr.Aeq=Aeq;
    dbgstr.beq=beq;
    dbgstr.Tgrid=Tgrid;
end


end

function [f,gradf] = objJ_fungrad(x, Jobj)
f = Jobj.J_Eval(x);
% Gradient of the objective function:
if nargout  > 1
    gradf = Jobj.GradJ(x);
end
end

function [f,gradf] = objft_fungrad(x, fobj)
f = fobj.feval(x);
% gradf = gradf_(fobj.ntm1+1:end);
% Gradient of the objective function:
if nargout  > 1
    gradf = fobj.fgrad(x);    %     hesf = hesf_(fobj.ntm1+1:end,fobj.ntm1+1:end);
end

end


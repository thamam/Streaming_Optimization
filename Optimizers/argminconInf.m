function [xstar] = argminconInf(farr,x0,MAXBUFSIZE, Ts)
%UNTITLED3 Summary of this function goes here
% JTinfObj{1,1} = T ; % 'Saved space - ptr to last row with with data% b -
% buff size JTinfObj{1,2} = bksz; x = (x_Tmb | x_Tmbp1,...,x_T); x_Tmb is
% not valid for update and provide as a refrence for the optimization
% proram Example: for T = 6, J= f2(x1,x2,x3,x4)+...+f6(x9,x10,x11,x12) x =
% [x1, x2 | x3, x4, ..., x9, 10| x11, x12] , x1-10 are the sol from T-1,
% and such that x1,x2 are frozen variables that are provided for the equality
% constraints of f2. x11,x12 are some arbitrary input variables to use for
% the intiialization, usually set as [1,1].

%Initialize 
A=[]; b=[]; Aeq=[]; beq=[];
bksz = farr{1,2};
T = farr{1,1};
Tmbp1 = max(2,T-MAXBUFSIZE+1);

% Determine if trucation is active
trncEn = (T>MAXBUFSIZE);

if T==2
    y0=x0; %for the first round only - no pre-solver
    
    f2 = farr{2,1};
    %prepare inequality constraints
    xInd = f2.ltObj.xtm1_ind(1):f2.ltObj.xt_ind(2);
    supp = f2.ltObj.f_sup; % [f2.ltObj.Tlim(1), f2.ltObj.Tlim(2)];
    Tgrid = supp(1):Ts:supp(2);
    A = -f2.ltObj.step(Tgrid,xInd).';
    b = zeros(size(A,1),1);
    
    [xstar] = fmincon(@(x) f2.feval(x),y0,A,b);
    
    
else % T >=3
    
    %% compute y0_T
    % we use pre-solver to compute y0_T=argmin_w(f_T(u,v)) s.t.
    % u=xhat_Tm1|Tm1 prepare fT for opt program
    fT =    farr{T,1};
    xInd_T = [fT.ltObj.xtm1_ind,fT.ltObj.xt_ind];
    supp_T = fT.ltObj.f_sup;
    Tgrid_T = supp_T(1):Ts:supp_T(2);
    A_T = -fT.ltObj.step(Tgrid_T,xInd_T).';
    b_T = zeros(size(A_T,1),1);
    y0_Tm1 = x0(end-2*bksz+1:end-bksz);
    Aeq_T = [eye(bksz),zeros(bksz)];
    beq_T = y0_Tm1;
    w0 = x0(end-2*bksz+1:end);
    y_tmp = fmincon(@(y) fT.feval(y),w0,A_T,b_T,Aeq_T,beq_T);
    y0=x0;
    y0(end-bksz+1:end) = y_tmp(end-bksz+1:end);
    
    %% compute y*
    %compute fun handle for last five functions in array
    Jarr = copyfromcell(farr,Tmbp1,T);
    N = min(MAXBUFSIZE,T-1);
    %Note that we use x as input but we force the first block to stay unchanged
    %     {1:min(MAXBUFSIZE,T-1)} = cell(farr{Tmbp1:T,1});
    
    
    %% prepare equality constraints
    if trncEn %(only when trucation is active)       
        Aeq = [eye(bksz), zeros(bksz, bksz*(min(T,MAXBUFSIZE+1)-1))];
        beq = y0(1:bksz);
    end
    
    %prepare inequality constraints
    xInd = Jarr{1}.ltObj.xtm1_ind(1):Jarr{end}.ltObj.xt_ind(2);
    supp = [Jarr{1}.ltObj.f_sup(1), Jarr{end}.ltObj.f_sup(2)];
    Tgrid = supp(1):Ts:supp(2);
    A = -Jarr{1}.ltObj.step(Tgrid,xInd).';
    b = zeros(size(A,1),1);
    
    % Call for fmincon solver
    [xstar] = fmincon( @(x) f_Eval(x,Jarr,N,bksz),y0,A,b,Aeq,beq);
end
end


function f = f_Eval(x,Jarr, N, bksz)
NT=numel(x);
f= 0; %@(x) 0;
for t=1:N
    Ut = computelocalproj(t,bksz,NT);
    f  = f + Jarr{t}.feval(Ut*x);
end

end

function [Ut] = computelocalproj(t,bksz,NT)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
U = eye(NT);
%             UiInd = [obj.fObjArray.datalist{t}.xtm1_idmap, ...
%                 obj.fObjArray.datalist{t}.xt_idmap];
UiInd = (t-1)*bksz+(1:2*bksz);
Ut = U( UiInd , : );
end
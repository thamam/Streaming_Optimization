function [step, l,u] = ONA_FBS(f0,H, bksz)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
l=0; 
u=0; 
step=0;
N = numel(f0);
d = N/bksz;



%% Backwards Solve  L U s_c =  -F(X_c)
Q{1} = H(1:bksz,1:bksz) ;  g{1} = -f0(1:bksz);
v{1} =   Q{1}\g{1};  %$ v_0 := Q_0^{-1} c_0 $
for  t=2:d     
    [g{t}, E{t-1}, H_t ] = retrieveBlocks(f0 ,H , bksz, N, t);
    
    U{t-1} =Q{t-1}\(E{t-1}.') ;    % 
    Q{t}  = H_t -  E{t-1} * U{t-1}; %
    
    v{t} = (Q{t})\ (g{t} - E{t-1}*v{t-1});   % 			 solve : $ v_k = Q_k^{-1}\left(c_k - E_k v_{k-1} \right) $
end
%% SOlve  U s = v  {Backward substitution}- \\
% 		All we need to do now, is to fold backward while solving the system $ Us = v $, where the Upper triangular matrix, U, and the intermediate solution vector $ v $, were obtained during the previous stage. 
step = v{d};     %   \s_{t_f} := v_{t_f} $
stp1 = step;
for t = (d-1):-1:1 
    st = v{t} - U{t}*stp1; % 
    stp1 = st;
    step = [st;step];
end

if nargout>1
%     U = eye(N);
    
    u=eye(N);
    l = zeros(N);
    l(1:bksz,1:bksz) = Q{1};    
    %construct lu
    for t=2:d
        bksft = (t-1)*bksz;
        UiInd = bksft +(1:bksz);
        
        %update upper triangular matrix
        u(UiInd-bksz, UiInd ) = U{t-1};
        
        %update lower triangular matrix
        l(UiInd,UiInd) = Q{t};
        l(UiInd,UiInd-bksz) = E{t-1};
    end
end

end


function [g_t, E_tm1,H_t] = retrieveBlocks(f0 ,H , bksz, N, t)
    U = eye(N);
    bksft = (t-1)*bksz;
    UiInd = bksft +(1:bksz);
%     Ut = U( UiInd , : );
    
    g_t = -f0(UiInd);
    H_t = H(UiInd,UiInd);
    E_tm1 = H(UiInd,UiInd-bksz);
 
end

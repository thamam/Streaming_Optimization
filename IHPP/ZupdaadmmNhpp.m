function theta = ZupdaadmmNhpp(xk,yk,N,P,gmap,rho  )
X = zeros(N,P);
%indicator matrix on the subset local variables x_i , W(i,j)=1 if (x_i)_j
%is active and zero o.w.
W=X;
for i=1:P
    X(gmap{1,i},i) = xk{1,i} +yk{1,i}/rho;
    W(gmap{1,i},i) = 1;
end

wv = sum(W,2);
theta = sum(X,2)./wv;
end
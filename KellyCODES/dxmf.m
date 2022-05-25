function dxu = dxmf(u)
%
% matrix-free partial derivative wrt x
% homogeneous Dirichlet BC
%
n2=length(u);
n=sqrt(n2);
h=.5*(n+1);
%
% turn u into a 2D array with the BCs built in
%
uu=zeros(n+2,n+2);
vv=zeros(n,n);
vv(:)=u;
uu(2:n+1,2:n+1)=vv;
%
% compute the partial
%
dxuu=zeros(n,n);
dxuu=uu(3:n+2,2:n+1)-uu(1:n,2:n+1);
%
% fix it up
%
dxu=h*dxuu(:);
end

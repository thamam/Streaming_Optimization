function dyu = dymf(u)
%
% matrix-free partial derivative wrt y
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
dyuu=zeros(n,n);
dyuu=uu(2:n+1,3:n+2)-uu(2:n+1,1:n);
%
% fix it up
%
dyu=h*dyuu(:);
end

grpcols = ceil(N/P);
% firstgrp = [1,2:1:m*(grpcols-1)+1];
subdivision = zeros(1,N);
subdivision(1)=1;
startind  = [2,2+m:m*(grpcols):(N-1)];
for i=1:numel(startind)
    curgrpind = startind(i)+(0:1:grpcols*m-1);
    subdivision(curgrpind)=i;
end



function [ map,weights,Kg ] = distmapping( options,C )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% We assign columns of C such that  each agent is assigned over entire
% colomun or columns. Since the first and the last column are the
% terminals (one node), the first and last agent will in general have
% different partition size assigned to them.

N=options.N;  P=options.P;  n=options.n;  m=options.m;

grpcols = ceil(n/P);
firstgrp = [1,2:1:m*(grpcols-1)+1];
subdivision = zeros(1,N);
subdivision(firstgrp)=1;
startind  = 2+m:m*(grpcols):N;
for n=1:numel(startind)-1
    curgrpind = startind(n)+(0:1:grpcols*m-1);
    subdivision(curgrpind)=n+1;
end

lastgrp = startind(end):1:N;
subdivision(lastgrp) = P;
map = repmat(subdivision,[N,1]);
weights = zeros(N,P);
for p=1:P
    CTC = C+C.';
    CTCi = reshape(CTC(map==p),N,[]);
    weights(:,p) = logical(sum(CTCi,2));
end

Kg=sum(weights,2);
% c = kron((1:options.P),ones(options.N,d));
% map = c(1:options.N,1:options.N);

end


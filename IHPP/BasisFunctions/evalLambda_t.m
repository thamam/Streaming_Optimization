function [ lambda_t ] = evalLambda_t( funcStruct,t,alpha )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% funcStruct.Type = 'haar';
J = funcStruct.J;
N = funcStruct.N;
hfunc = funcStruct.fHandle;

phi_n_time = zeros(numel(t),N);
for n=0:(N-1)
    phi_n_time(:,n+1) = hfunc(t,n,J).';
end

lambda_t = phi_n_time * alpha;


end


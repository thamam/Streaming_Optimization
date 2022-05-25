%quadpenfromxbar
clear;
close all;
clc;
%% Description
% We compute in addition to the constrained also unconstrained with penalty
% term given by the distance to the feasible set
%% Define  the two objective quadratics
N=2;
[Q1, egsQ1] = generateQuadratic(N);
[Q2, egsQ2] = generateQuadratic(N);


% Lift and join the matrices
EYE = eye(3);
U1 = EYE(1:3,1:2);
U2 = EYE(1:3,2:3);
F1 = U1*Q1*U1.';
%%
F2 = U2*Q2*U2.';

F= F1+F2;
[Uf,Lambda] = eigs(F)
x3 = Uf(:,3);

u3 = x3(1:2);
v3 = x3(2:3);

a = u3.'*Q1*u3;
b = v3.'*Q2*v3;

u3nrm = u3.'*u3;
v3nrm = v3.'*v3;

a+b

Eq = [min(eigs(Q1)), min(eigs(Q2))]
sum(Eq)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===============================================================
function [Xc] = proj2C(x, A, b)
    Xc =  fmincon(@(y) norm(x-y)^2, x,A,b);
end

function [Q, Qeis] = generateQuadratic(N)
Q_ = randn(N,N);
Q_ = Q_*Q_.';
[V, ~] = eig(Q_*Q_.');
d = abs(0.1*randn(N,1));
Qu = V*diag(d)*V.';
eQ = eigs(Qu);
meQ = sum(eQ)/N;
Q = Qu/meQ;
Qeis = eigs(Q);
end

function [fout] = foutarray(f,X,Y)
    fout = zeros(size(X));
    [r,c]=size(X);
    for i=1:r
        for j=1:c
            fout(i,j) = f([X(i,j);Y(i,j)]);
        end
    end           
end
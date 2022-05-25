%  Level 1
function [xstar,outputArg2] = OFSol(inputArg1,inputArg2)
%OPTFILTSOLVER Summary of this function goes here
%   Detailed explanation goes here

isBar=1; alpha0=1 ; mu = 10; tol;

K =  [];%number of constraints
if isBar %Solving with log barrier method
    alpha0 = [];
    y=y0;
    while K/alpha >= tol
        alpha = alpha0/mu;        
        lgbrmodify(fObj, alpha)% Modify objective arguments       
        % solve
        [xstar] = FBNalg(y,fObj,algParams)
        y=xstar;
    end
else %Direct solver
    [xstar] = FBNalg(y,fObj,algParams)

    
end


end


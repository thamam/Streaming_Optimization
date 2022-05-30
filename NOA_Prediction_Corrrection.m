function [varargout] = NOA_Prediction_Corrrection(nhp_params,streamOptions,tK,bksz)
% NOA
% [xstrm, XHIST, streambrkdwnData] = NOA(nhp_params,streamOptions, tK)
%  xstrm -  as the optimal streaming solution including also the
%       variables that were fixed after max buffer size was reached
%  XHIST - Solution array where each column corresponds to the "active"
%       estimates at time T=1,2,.... If variables exceeded memory size, they
%       will have the value 0
%  streambrkdwnData - contains breaking down tables of the streaming
%       batches
%% Truncated Online Newton Algorithm (NOA)

%% Extract Input parmeters


lb_delta=streamOptions.lb_delta;
nonegxMod =streamOptions.nonegxMod; 
nhp_params.blocksize = bksz;
nhp_params.funcNum = nhp_params.d/bksz-1;
nhp_params.blocksNum = nhp_params.funcNum + 1;
nhp_params.fminconOptions = ...
    optimoptions(@fmincon,'MaxIterations',300,'MaxFunctionEvaluations',50000, ...
    'OptimalityTolerance',1e-10,'ConstraintTolerance',1e-10);
fprintf(' $$ Setting streaming window size \n\n');
[Tlim_div, xidmapTbl , xIndTbl, xlimTbl] = streamingbreakdown(nhp_params);  % Screen printout
streambrkdwnData  = struct('Tlim_div',Tlim_div, 'xidmapTbl',xidmapTbl,...
    'xIndTbl',xIndTbl, 'xlimTb',xlimTbl);
Tlim_div=[0,0;Tlim_div]; % first row is zeros because T starts from 1 and blocks start from 0
Tend = size(Tlim_div,1);

%% Initialization
% JTObj     = JT_BuffObjCLASS(MAXBUFSIZE, bksz);
% JTObj= JT_BuffObjCLASS(MAXBUFSIZE, bksz);
% xbufObj   = MEMBUFFERCLASS(MAXBUFSIZE+1,bksz,'double',true); % buffer to manage active xvariables - \xtld_t
% size set to MAXBUFSIZE+1 to match function buffer which has a size of
% MAXBUFSIZE

% JTinfObj  = cell(Tend,1);
% do not change data written to frist row of array - saved for info on
% array

% JTinfObj{1,2} = bksz;
x0=(1:1:bksz).'*0+1;

%% Temporary location before refactoriing
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\
%% Prediction-Correction code
Xpc=[];

%% Initialize - during first solve, g_1 = f(x0, x1) or get initial solver
isfirstround = 1;
% Provide initial guess for (x0,x1)
xT0=ones(1,1);

for  T=2:Tend
    xidmap = {xidmapTbl(T-1,:), xidmapTbl(T,:)}; % Opt. variable ids
    xInd ={xIndTbl(T-1,:), xIndTbl(T,:)};       % the spline indices corresponding to the opt. variable ids
    nhpObj_T = NHP_ftObjCLASS(nhp_params.splineOrder, tK, Tlim_div(T,:), xInd,lb_delta);
    fT = ftObjCLASS(nhpObj_T,T,xidmap);    % Input - new increment of f_T(x_{T-1},x_T)         %

    %       OPTIONAL for Tt  % use finer stepping inside Tlim_div(T,:) window
    %         % define current cost function
    if isfirstround
        x0 = ones(2*bksz,1);
        gt= @(x) obj_pc_fungrad(x,fT, true);
    else
%         x0 = xT0;
%         gt= @(x) obj_pc_fungrad(x,fT, false, Xpc(end));
        gt= @(x) obj_pc_fungrad(x,fT, true);
        x0 = [Xpc(end); xT0];

    end

    % correction step
    isham = 1 ;  rsham = 0 ; % is Newton's method
    maxit = 300; % maxmium number of nonlinear iterations
    atol = 1.d-6; rtol = 1.d-6;
    t0 = lb_delta; mu=10; eps=1e-6;
    tol = [atol, rtol] ; %relative/absolute tolerance
    parms = [maxit, isham, rsham];
    options.UseHessian=1;
    options.bksz = bksz;
    %         [xhat_cor,iterHist] = NwBarrier(x_pred, gt, parms, options, tol,t0, mu, eps);
    [yout, it_hist] = NwBarrier(x0, gt, parms, options, tol, t0, mu, eps);

    Xpc([xidmapTbl(T-1,:), xidmapTbl(T,:)],T) = yout;

%     if isfirstround
%     else
%         Xpc = [Xpc; yout(end)];
%     end
% 
%     if isfirstround
%         isfirstround = 0;
%     end
end
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
varargout{1} = Xpc;
if nargout>1
    varargout{2} = it_hist;
end

end


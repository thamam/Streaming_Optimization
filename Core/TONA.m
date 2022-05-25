function [varargout] = TONA(nhp_params,streamOptions,tK, MAXBUFSIZE)

% [xstrm, XHIST, streambrkdwnData] = TONA(nhp_params,streamOptions, tK)
%  xstrm -  as the optimal streaming solution including also the
%       variables that were fixed after max buffer size was reached
%  XHIST - Solution array where each column corresponds to the "active"
%       estimates at time T=1,2,.... If variables exceeded memory size, they
%       will have the value 0
%  streambrkdwnData - contains breaking down tables of the streaming
%       batches
%% Truncated Online Newton Algorithm (TONA)

%% Extract Input parmeters

% MAXBUFSIZE = streamOptions.MAXBUFSIZE;
lb_delta=streamOptions.lb_delta;
nonegxMod =streamOptions.nonegxMod;
bksz = streamOptions.bksz;


%%


nhp_params.blocksize = bksz;
nhp_params.funcNum = nhp_params.d/bksz-1;
nhp_params.blocksNum = nhp_params.funcNum + 1;
nhp_params.fminconOptions = ...
    optimoptions(@fmincon,'MaxIterations',300,'MaxFunctionEvaluations',50000, ...
    'OptimalityTolerance',1e-10,'ConstraintTolerance',1e-10);

[Tlim_div, xidmapTbl , xIndTbl, xlimTbl] = streamingbreakdown(nhp_params);
streambrkdwnData  = struct('Tlim_div',Tlim_div, 'xidmapTbl',xidmapTbl,...
    'xIndTbl',xIndTbl, 'xlimTb',xlimTbl);
Tlim_div=[0,0;Tlim_div]; % first row is zeros because T starts from 1 and blocks start from 0
Tend = size(Tlim_div,1);

%% Initialization
JTObj     = JT_BuffObjCLASS(MAXBUFSIZE, bksz);
% JTObj= JT_BuffObjCLASS(MAXBUFSIZE, bksz);
xbufObj   = MEMBUFFERCLASS(MAXBUFSIZE+1,bksz,'double',true); % buffer to manage active xvariables - \xtld_t
% size set to MAXBUFSIZE+1 to match function buffer which has a size of
% MAXBUFSIZE

JTinfObj  = cell(Tend,1);
% do not change data written to frist row of arracy - saved for info on
% array 

JTinfObj{1,2} = bksz;
x0=(1:1:bksz).'*0+1;
xhat_strm=x0;
Xflushed=[]; %xflushed records all the optimization variable that were discarded from the active memory(buffer)
%% Main loop - simulating streaming arrival of batches
tic
for T=2:Tend
    JTinfObj{1,1} = T ; % 'Saved space - ptr to last row with with data
    %==========================================================================================%
    % Instantiate fT and add it to the array
    % %
    %==========================================================================================%
    xidmap = {xidmapTbl(T-1,:), xidmapTbl(T,:)};                                               %
    xInd ={xIndTbl(T-1,:), xIndTbl(T,:)};                                                      %
    nhpObj_T = NHP_ftObjCLASS(nhp_params.splineOrder, tK, Tlim_div(T,:), xInd,lb_delta);
    fT = ftObjCLASS(nhpObj_T,T,xidmap);    % Input - new increment of f_T(x_{T-1},x_T)         %
    JTObj.addfunction(fT,T);
    %save f_T to array of all f_t instantiations
    JTinfObj{T,1} = fT;
    JTinfObj{T,2} = T;
    %==========================================================================================%
    xT0=ones(fT.nt,1);
    x0=[xhat_strm;xT0];
    
    %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    %% Run streamlined solver
    options.nonegxMod = nonegxMod;
    options.Solver = 'NW'; 
    mu=10;
    eps=1e-6;
    %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    [xhat_NW , iterHist] =  NWSolver(JTObj,x0,lb_delta,mu, eps) ;
    xhat_strm = xhat_NW;
    %% Save xhat to hist array
    indShift = bksz*max(0,T-(MAXBUFSIZE));
    maxind = min((T),MAXBUFSIZE);
    histInd = indShift + (1:bksz*maxind);
    XHIST(histInd,T) = xhat_strm;
    
    %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\
    % Option for batch run - [0, Tcurrent]
    % TBD - take out of this code and implement as a seperate code
%     if false 
%         
%         [xbatch, fval_batch, DBG_batch] = batchSolver(Tlim_div,tK, xIndTbl,...
%             T, xidmapTbl, nhpObj_T, spord, Ts, nhp_params, options);                
%     end
    
    %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    if T>=MAXBUFSIZE
        Xflushed=[Xflushed;xhat_strm(1:bksz)];
    end    
end

xstrm = XHIST(:,end);
if MAXBUFSIZE<numel(xstrm) %test that truncation was performed
    xstrm(1:numel(Xflushed))=Xflushed; %retrieve \xf_t from memory
end

varargout{1} = xstrm;
if nargout >1
    varargout{2} = XHIST;
    if nargout>2
        varargout{3} = streambrkdwnData;
    end
end


end


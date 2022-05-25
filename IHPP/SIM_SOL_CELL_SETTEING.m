%% SIM_SOL_CELL_SETTEING
s=0;
%% Gloal Solution
sim  = struct  ;
sim.Type = 'GlobalCvx';
SIM{s+1}=sim;
s=s+1;
sectionLength = params.sectionLength;


%% Waves Propogation
if 0
    sim  = struct  ;
    sim.Type = 'wavesCD';
    WCD_opts = struct ('sectionLength',sectionLength,'MaxIter',100,'ActivationMode','WavesProp' );
    sim.options = WCD_opts;
    SIM{s+1}=sim;
    s=s+1;
end

if false
    %% FWD_W_single_BWD_UPDATE(JR)
    % For each new frames section: 1. solve new frame fixing edge coeffs
    % 2. propgate solution backword only :
    % \tau_F = -1,
    % \tau_R = 1
    
    sim  = struct  ;
    sim.Type = 'FWD_STRM_BWD_UPDATE_JR';
    WCD_opts = struct ('sectionLength',sectionLength,'MaxIter',[],'ActivationMode','FWD_STRM_BWD_UPDATE_JR' );
    W= ceil(params.duration/sectionLength);  % number of sections
    WCD_opts.MaxIter = 2*W-1;
    sim.options = WCD_opts;
    SIM{s+1}=sim;
    s=s+1;
    
end
%% Coordinate Descent
if false
    sim  = struct  ;
    sim.Type = 'CoordinateDescent';
    WCD_opts = struct ('sectionLength',sectionLength,'MaxIter',500,'ActivationMode','parallel_CD');
    sim.options = WCD_opts;
    SIM{s+1}=sim;
    s=s+1;
end

%% Sweeping FWD and BWD
if false
    % \gamma_G = 1
    % \gamma_L = 1
    % \tau_F=tau_R = 1;
    %gamma_F = gamma_R = 0;
    sim  = struct  ;
    sim.Type = 'SWEEP';
    WCD_opts = struct ('sectionLength',sectionLength,'MaxIter',500,'ActivationMode','SWEEP');
    sim.options = WCD_opts;
    SIM{s+1}=sim;
    s=s+1;
end
%% Random coordinate descent
% At each iteration pick at random which coordinate to solve
% % sim  = struct  ;
% % sim.Type = 'RANDOMCD';
% % WCD_opts = struct ('sectionLength',sectionLength,'MaxIter',500,'ActivationMode','RANDOMCD');
% % sim.options = WCD_opts;
% % SIM{s+1}=sim;
% % s=s+1;

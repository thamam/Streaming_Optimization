% Selected Results for display

%% Coordinate descent 
% lambda is 45, time from [2,69], each section is 4s
clear ;clc
close all;
load('results_W-06_08_16_T- 4_52_PM_GCD_simData_maxRate_45_T_2_69_L_3.mat');
saveVideo= true;
plotScript_coordniateDescent

%%
shtdw = 'n';% input('Do you want to restart your computer? \n','s')
if strcmpi(shtdw,'y')
     command = 'shutdown -r -t 0';
    [status,cmdout] = system(command); 
end

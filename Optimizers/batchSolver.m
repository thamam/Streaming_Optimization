function [xbatch, fval_batch, DBG_batch] = batchSolver(Tlim_div,tK, xIndTbl, T, xidmapTbl,...
    nhpObj_T, spord, Ts, nhp_params, options)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



Tlim_batch = [Tlim_div(2,1),Tlim_div(T,2)];
tK_batch = tK(tK >= Tlim_batch(1) & tK <= Tlim_batch(2));
Tmed= round(numel(T)/2);
xidmapBatch = { xidmapTbl(1):1:xidmapTbl(Tmed,2),...
    xidmapTbl(Tmed+1,1):1:xidmapTbl(T,2)};
xIndBatch_carry ={ xIndTbl(1):1:xIndTbl(Tmed,2),...
    xIndTbl(Tmed+1,1):1:xIndTbl(T,2)};
xInd_batch = cell2mat(xIndBatch_carry);
support_batch = minmax(reshape(nhpObj_T.support([xInd_batch(1);xInd_batch(end)]),1,[]));
nhpObj_batch_T = NHP_ftObjCLASS(spord, tK_batch, Tlim_batch, xIndBatch_carry);                                 %
fbatch_T = ftObjCLASS(nhpObj_batch_T,T,xidmapBatch);
x0fs = xInd_batch.'*0+1;
%
[xbatch, fval_batch, DBG_batch] = argminBatch(x0fs,fbatch_T,...
    cell2mat(xIndBatch_carry),support_batch,Ts, nhp_params.fminconOptions , options.nonegxMod);%, Tlim_batch );
%%
%         [xNw_batch, Fbatch, Hbatch] = NwBatchSolver(x0fs,fbatch_T,...
%             cell2mat(xIndBatch_carry),support_batch,Ts, options,lb_delta, mu, eps);%, Tlim_batch );



end


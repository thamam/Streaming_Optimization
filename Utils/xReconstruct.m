function [XREC] = xReconstruct(XHIST,bksz,MAXBUFSIZE)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
Tend = size(XHIST,2);
XREC = XHIST(:,1:1+MAXBUFSIZE);
lastInd = 0 ;
for t=MAXBUFSIZE+bksz:Tend
    lastInd = lastInd + bksz;
    XREC(:,t) = XHIST(:,t);
    XREC(1:lastInd,t) = XREC(1:lastInd,t-1);
end


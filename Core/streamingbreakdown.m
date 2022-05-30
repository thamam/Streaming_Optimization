function [Tlim_div, xidmap , xindTbl, xindlimTbl] = streamingbreakdown(params)
%streamingbreakdown Summary of this function goes here
%   Detailed explanation goes here
blocksize = params.blocksize;
Tlim = [params.simTstart, params.simTfinish];
%
%%% Create division to time blocks
%
% L = (params.splineOrder+1)/2;
%
%%% Create division of basis function endpoints
%
startIndCol= (params.nInd(1):blocksize:params.nInd(end)).';
xindlimTbl = [startIndCol , startIndCol + (blocksize-1)];

numofblks=ceil(params.d/blocksize);
numoffuncts = numofblks - 1;
xindTbl=zeros(numofblks,blocksize) ;
%Expand the nInd in each local nlim
for i=1:numofblks
    xindTbl(i,:) = xindlimTbl(i,1):1:xindlimTbl(i,2);
end
%compute the trip point - finding the transition point from block
%[x_t-1,x_t] to [x_t,x_{t+1}];
trippoint=zeros(numoffuncts-1,1);
for i=2:numoffuncts
    [trippoint(i-1), ~] =  compuBasisSupp(Tlim, params.splineOrder, xindlimTbl(i+1,1));
end


Tlim_div = [[Tlim(1);trippoint(:)], [trippoint(:);Tlim(2)]];
%
%%% Expand nind within each time frame
%
%
%%% define mapping for x from local vector to global vector
xidmap=reshape(1:params.d, [], size( xindTbl,1)).';
% xmap=reshape(1:numel(nindDiv), [], size( nindDiv,1)).';
end


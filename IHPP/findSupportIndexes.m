function [ nInd ] = findSupportIndexes( phi_n,start,finish,L )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% fimd minimal index for which phi_n is supported
    done=[0,0];
    lowInd = ceil(start-(L+1)/2);
    upInd = floor(finish+(L+1)/2);
    while(done(1)~=1)
        if(phi_n(start,lowInd)==0)
            lowInd = lowInd+1;
        elseif (phi_n(start,lowInd)~=0 && phi_n(start,lowInd-1)~=0)
            lowInd = lowInd-1;
        else
            done(1)=1;
        end         
    end
    
    while(done(2)~=1)
        if(phi_n(finish,upInd)==0)
            upInd = upInd-1;
        elseif (phi_n(finish,upInd)~=0 && phi_n(finish,upInd+1)~=0)
            upInd = upInd+1;
        else
            done(2)=1;
        end         
    end
    
    nInd=max(lowInd,floor(start)):1:min(ceil(finish), upInd);


end


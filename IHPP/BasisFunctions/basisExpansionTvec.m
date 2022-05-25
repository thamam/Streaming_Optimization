function [ f_t ] = basisExpansionTvec( hfunc,t )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    for i=1:numel(t)
        f_t(i) = hfunc(t(i));
    end
 
end


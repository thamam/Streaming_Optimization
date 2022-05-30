function [Jobj]=updatearrayprop(Jobj, pval, Tmbp1, MAXBUFSIZE)
%update property value over a cellarray of class
% obj is handle to the class object cell array
% pname is the name of the property field to be upodated
% pval is the new value to be assigned
for i=Tmbp1:MAXBUFSIZE
    propmodcom = sprintf('Jobj.fObjArray.datalist{%i}.ltObj.lb_delta=pval ;',i);
    eval(propmodcom);
end
end
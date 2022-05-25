function [ indMatList ] = extractIndList( supMat )
% Returns explicit indices for alpha indices included within the support of
% each frame
[W,~]=size(supMat);
indMatList = [];
for w=1:W
    indMatList(w,:) = supMat(w,1):supMat(w,2);
end
end


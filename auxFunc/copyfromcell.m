function [Src] = copyfromcell(Tgt,istart,iend)
    

    cnt=0;
    for k=istart:iend
        cnt=cnt+1;
        Src{cnt} = Tgt{k};
    end

end
classdef JTObjCLASS<handle
    %FOBJCLASS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nmfnc %total expected number of batches for memory allocation
        bfsz % buffer size
        bksz % block size
        fObjArray % cell array of all the currently sotred ft objects
        T
        NT %total number of variables = sum_t nt
        Tlim, xIndtbl, tKCnt
        
        
    end
    
    methods
        function obj = JTObjCLASS(nmfnc, bfsz, bksz)
            %JTObjCLASS Construct an instance of this class
            %   Detailed explanation goes here
            obj.nmfnc = nmfnc ;
            obj.bfsz = bfsz;
            obj.bksz = bksz;
            obj.fObjArray{1}=0;
            obj.NT = bksz;
            obj.Tlim=[0,0];
            obj.tKCnt = [];
                        
        end
        
        function obj = addfunction(obj,fT,T)
            %METHOD1 Add an objective
            obj.fObjArray{T} = fT;
            obj.T=T;
            obj.NT = obj.NT+fT.nt;
            if T==2
                obj.Tlim(1) = fT.ltObj.Tlim(1);
                obj.xIndtbl(T-1,:)=fT.ltObj.xtm1_ind;
            end
            obj.Tlim(2) = fT.ltObj.Tlim(2);
            obj.xIndtbl(T,:)=fT.ltObj.xt_ind;
%             obj.tKCnt = [obj.tKCnt; fT.ltObj.tK_t];
            obj.tKCnt = obj.tKCnt + numel(fT.ltObj.tK_t);
        end
        
        function J = J_Eval(obj,x)
            NT=numel(x);
            J= 0; %@(x) 0;
            for t=2:obj.T
                Ut = computelocalproj(obj, t);
                J  = J + obj.fObjArray{t}.feval(Ut*x);
            end
        end
        
        function [Ut] = computelocalproj(obj, t)
            %UNTITLED2 Summary of this function goes here
            %   Detailed explanation goes here
            
            U = eye(obj.NT);
            
            UiInd = [obj.fObjArray{t}.xtm1_idmap,obj.fObjArray{t}.xt_idmap];
            Ut = U( UiInd , : );
            
            
        end
        
    end
end


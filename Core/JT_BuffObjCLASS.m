classdef JT_BuffObjCLASS<handle
    %JT_BuffObjCLASS Comulative objecive function object: 
    %   This object aggregates all the ft funcational objects as they are
    %   obtained. It's main methods are evaluating the cost function, its
    %   gradient and its Hessian
    
    properties
        nmfnc %total expected number of batches for memory allocation
        bfsz % buffer size
        bksz % block size
        fObjArray % cell array of all the currently sotred ft objects
        T %counter index of received bataches 
        NT % Nt=\sum_{t=1}^T R^{n_T}
        Tlim, xIndtbl, tKCnt
        support % the maximal domain for which xInd(i) \neq 0 for some i and at some t
    end
    
    methods
        function obj = JT_BuffObjCLASS(bfsz, bksz)
            %JT_BuffObjCLASS Construct an instance of this class
            %   Detailed explanation goes here
            obj.bfsz = bfsz;
            obj.bksz = bksz;
            obj.fObjArray = MEMBUFFERCLASS(bfsz,1,'ftObjCLASS');
            obj.NT = bksz;
            obj.Tlim=[0,0];
        end
        
        function obj = addfunction(obj,fT,T)
            %METHOD1 Add an objective
            fTmb = obj.fObjArray.add(fT);
            obj.T=T;
            obj.NT = min(obj.NT+fT.nt,(obj.bfsz+1)*obj.bksz); %bfs+1 since numblkc=num_of_funct+1
            if T==2
                obj.Tlim(1) = fT.ltObj.Tlim(1);
                obj.xIndtbl(T-1,:)=fT.ltObj.xtm1_ind;
                obj.xIndtbl(T,:)=fT.ltObj.xt_ind;
                obj.support = fT.ltObj.f_sup;
            elseif T>2 && T<=obj.bfsz+1 % while buffer truncating isn't active
                obj.xIndtbl(T,:)=fT.ltObj.xt_ind;
                obj.support(2) = fT.ltObj.f_sup(2);
            else % buffer is truncating and buffer info need both add and subtraction of data info
                obj.Tlim(1) = obj.fObjArray.datalist{1}.ltObj.f_sup(1);
                xIndTbltmp = [obj.xIndtbl; ...
                    fT.ltObj.xt_ind];
                obj.xIndtbl = xIndTbltmp(2:end,:);
                obj.support(1) = obj.fObjArray.datalist{1}.ltObj.f_sup(1);
                obj.support(2) = fT.ltObj.f_sup(2);
            end
            % things that we do for every T
            obj.Tlim(2) = fT.ltObj.Tlim(2);
        end
        
        function J = J_Eval(obj,x)
            J= 0; 
            lptr = obj.fObjArray.lptr;
            N = min(obj.T-1,obj.bfsz);
            for t=1:N
                Ut = computelocalproj(obj,t);
                % Ut is (ntm1+nt) times NT matrix with ones blocks in the
                % tm1 block row by first block columns and t row block  by
                % the second the block column.
                
                J  = J + obj.fObjArray.datalist{lptr-1+t}.feval(Ut*x);
            end
        end
        
        function gradJ = fgrad(obj,x)
            NT =  obj.NT; % Nt=\sum_{t=1}^T R^{n_T}
            gradJ= zeros(NT,1);
            lptr = obj.fObjArray.lptr;
            N = min(obj.T-1,obj.bfsz); % # of functions in buffer
            for t=1:N
                Ut = computelocalproj(obj,t); 
                % Ut is (ntm1+nt) times NT matrix with ones blocks in the
                % tm1 block row by first block columns and t row block  by
                % the second the block column.
                
                fgrad_t=obj.fObjArray.datalist{lptr-1+t}.fgrad(Ut*x);
                % lift fgrad_t to R^NT by zeros padding and combine W/
                % gradJ
                gradJ = gradJ + Ut.'*fgrad_t;                                
            end
        end
        
        function HessJ = fhessian(obj,x)
            NT =  obj.NT; % Nt=\sum_{t=1}^T R^{n_T}
            HessJ= zeros(NT, NT);
            lptr = obj.fObjArray.lptr;
            N = min(obj.T-1,obj.bfsz); % # of functions in buffer
            for t=1:N
                Ut = computelocalproj(obj,t); 
                % Ut is (ntm1+nt) times NT matrix with ones blocks in the
                % tm1 block row by first block columns and t row block  by
                % the second the block column.
                
                fHess_t=obj.fObjArray.datalist{lptr-1+t}.fhessian(Ut*x);
                % lift fgrad_t to R^NT by zeros padding and combine W/
                % gradJ
                HessJ = HessJ + Ut.'*fHess_t*Ut;                                
            end
        end
        
        function [Ut] = computelocalproj(obj, t)
            %UNTITLED2 Summary of this function goes here
            %   Detailed explanation goes here
            U = eye(obj.NT);
            %             UiInd = [obj.fObjArray.datalist{t}.xtm1_idmap, ...
            %                 obj.fObjArray.datalist{t}.xt_idmap];
            UiInd = (t-1)*obj.bksz+(1:2*obj.bksz);
            Ut = U( UiInd , : );
            
            
        end
    end
end


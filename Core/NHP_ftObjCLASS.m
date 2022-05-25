classdef NHP_ftObjCLASS < handle&nhpp_bspline
    %NHP_ftObjCLASS This is the class implementation for ft in the
    %nonhomogenous Poisson process. The vector form of ft()is given by: by
    % f(xtm1,xt) = <at,xt>+<bt,xtm1> - sum log( <c_(t,m),x_t> + <d_(t,m),xtm1>
    % )
    
    properties
        spord
        ntm1, nt %size of blocks xtm1 and xt
        xtm1_ind, xt_ind
        tK_t
        Tlim % limits of blocks I_t time interval
        f_sup % support of function
        at, bt, c_tm, d_tm
        
        %%% log barrier params %%%
        lb_delta  %log barrier weight
%         tCon  % interval grid where non-neg constraints imposed
        BtCon,AtCon, V_BA % (nxK matrix, with (At)_ij = phi_i(t_j) t_j in I_t
        K  % numel(tCon)
        lbgrMod =0;%is set to 1 if lgbr is used
        
    end
    
    methods(Access = public)
        function obj = NHP_ftObjCLASS(spord, tK, Tlim, xind,lb_delta)
            obj@nhpp_bspline(spord); %constructor of <nhpp_bspline
            obj.xtm1_ind = xind{1};
            obj.xt_ind = xind{2};
            tK_t = tK(tK>Tlim(1) & tK<Tlim(2));
            obj.tK_t = tK_t;
            obj.ntm1 = numel(obj.xtm1_ind);
            obj.nt = numel(obj.xt_ind);
            obj.Tlim=Tlim;
            obj.f_sup =[min(obj.support(obj.xtm1_ind(1))),  max(obj.support(obj.xt_ind(end)))];
            obj.spord=spord;
            %compute the terms that are independent of x
            at = zeros(obj.nt,1);
            bt = zeros(obj.ntm1,1);
            
            for i=1:obj.ntm1
                bt(i) = integral(@(t) obj.step(t,obj.xtm1_ind(i)),Tlim(1),Tlim(2));
            end
            
            for j=1:obj.nt
                at(j) = integral(@(t) obj.step(t,obj.xt_ind(j)),Tlim(1),Tlim(2));
            end
            
            c_tm = (obj.step(tK_t,obj.xt_ind(:)));
            d_tm = (obj.step(tK_t,obj.xtm1_ind(:)));
            
            obj.at= at;
            obj.bt=bt;
            obj.c_tm = c_tm;
            obj.d_tm = d_tm;
            
            %Log barrier initialization
%             if nargin>4
%                 obj.lbgrMod=1;
%                 obj.lb_delta = lb_delta;
%                 obj.Ts=Ts;
%                 tCon = (obj.f_sup(1)+Ts):Ts:(obj.f_sup(2)-Ts);
%                 obj.tCon = tCon;
%                 obj.K = numel(tCon);
%                 obj.BtCon = obj.step(tCon,obj.xtm1_ind(:));
%                 obj.AtCon = obj.step(tCon,obj.xt_ind(:));
%                 obj.V_BA = [obj.BtCon;obj.AtCon];
%             end
            
        end
        
        function fout = ft(obj, xtm1,xt)
            % fout returns f(xtm1,xt) = <at,xt>+<bt,xtm1> - sum log( <c_(t,m),x_t>
            % + <d_(t,m),xtm1> )
            cdTx = obj.c_tm.'*xt + obj.d_tm.'*xtm1;
            log_cd = log(cdTx);
            fobj = obj.at.'*xt + obj.bt.'*xtm1 - sum(log_cd); % Objective cost
            
%             if obj.lbgrMod
%                 lbOut = lbf_copm(obj,xtm1,xt); % log barrier cost
%                 fout =  obj.lb_delta*fobj + lbOut;
%             else
                fout =  fobj ;
%             end
        end
        
        function fgrad = grad(obj, xtm1,xt)
            % fgrad returns \nabla f(xtm1,xt) = at+bt-sum log(
            % \frac{c_(t,m)+d_(t,m)}{ <c_(t,m),x_t> + <d_(t,m),xtm1>} )
            
            aplusb = eye(obj.nt+obj.ntm1)*[obj.bt;obj.at];
            
            cplusd = eye(obj.nt+obj.ntm1)*[obj.d_tm;obj.c_tm];
            cdTx = (obj.c_tm.'*xt + obj.d_tm.'*xtm1).';
            
            fobjgrad = aplusb -  sum(cplusd ./ cdTx,2);
            
%             if obj.lbgrMod
%                 lb_grad_ = lbfgrad_copm(obj,xtm1, xt);
%                 fgrad = obj.lb_delta*fobjgrad + lb_grad_;
%             else
               fgrad = fobjgrad;
%             end
                
        end
         
        function [varargout] = hessian(obj, xtm1,xt )
            % Compute the Hessian Matrix
            
            cplusd = eye(obj.nt+obj.ntm1)*[obj.d_tm;obj.c_tm];
            cdTx = (obj.c_tm.'*xt + obj.d_tm.'*xtm1).';
            N = obj.ntm1+obj.nt;
            
            Pk_=cplusd.';
            den2 = (cdTx).^2 ; % den = [ sum(x.'*p0) , sum(x.'*p_1),...,sum(x.'*p_K)]
            Hout_=zeros(N);
            for i=1:N
                pi = Pk_(:,i).' ;
                for j=1:N
                    pj = Pk_(:,j).' ;
                    Hout_(i,j) = sum(pi.*pj./den2);
                end
            end
            
%             if obj.lbgrMod
%                 lb_Hes_f = lbfHessian_copm(obj,xtm1,xt, N);
%                 Hout = obj.lb_delta*Hout_ + lb_Hes_f;
%             else
                Hout = Hout_;
%             end
            
            %             if nargout ==2
            %                 [l,u] = lu(Hout);
            %                 varargout{1} = l;
            %                 varargout{2} = u;
            %             else
            varargout{1} = Hout;
            %             end
            
        end
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %  Log barrier functions
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function lb_f      = lbf_copm(obj,xtm1, xt)
            oneK_T=ones(1,obj.K);
            lb_f = -oneK_T*log([obj.BtCon.',obj.AtCon.']*[xtm1; xt]);
        end
        
        function lb_grad_f = lbfgrad_copm(obj,xtm1, xt)
            den = (obj.V_BA.'*[xtm1;xt]).' ; %den = [ sum(x.'*a_t0) , sum(x.'*a_t1),...,sum(x.'*a_tK)]
            lb_grad_f =   -sum(obj.V_BA ./ den,2);
        end
        
        function lb_Hes_f  = lbfHessian_copm(obj,xtm1, xt, N)
            V_BAT=obj.V_BA.';
            den2 = ((V_BAT*[xtm1;xt]).').^2 ; % den = [ sum(x.'*p0) , sum(x.'*p_1),...,sum(x.'*p_K)]
            Hout=zeros(N);
            for i=1:N
                BAt_i = V_BAT(:,i).' ;
                for j=1:N
                    BAt_j = V_BAT(:,j).' ;
                    Hout(i,j) = sum(BAt_i.*BAt_j./den2);
                end
            end
            lb_Hes_f = Hout;
        end
        
    end
    
end



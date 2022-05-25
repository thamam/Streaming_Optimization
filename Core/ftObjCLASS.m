classdef ftObjCLASS<handle
    %FTOBJCLASS This is the interface/wrapper for the basic loss term
    %object
    %--------  end of class help -------------%
    
    
    properties
        %         ft, Ft, Ht   % function and it's 1st and 2nd derivative handles
        t   % block id
        xtm1_idmap, xt_idmap % variable indices whithin the global vector \xbar_T
        Nt  % total size of xtm1 and xt blocks
        ltObj %the object of the implemented loss function
        nt
        ntm1
    end
    
    methods
        function obj = ftObjCLASS(ltObj,t,xidmap)
            %FTOBJCLASS Construct an instance of this class
            %   Detailed explanation goes here
            obj.ltObj = ltObj;
            obj.t = t;
            obj.xtm1_idmap = xidmap{1};
            obj.xt_idmap = xidmap{2};
            obj.Nt = numel(cell2mat(xidmap));
            obj.ntm1 = numel(xidmap{1});
            obj.nt = numel(xidmap{2});
        end
                
        function Ftag = fhessian(obj,varargin)
            
            if nargin==2
                x = varargin{1};
                midInd = obj.ntm1;
                xtm1 = x(1:midInd);
                xt = x(midInd+1:end);
            else
                xtm1 = varargin{1};
                xt = varargin{2};
            end
            
            Ftag =  obj.ltObj.hessian(xtm1,xt) ;
        end
                       
        function Ft = fgrad (obj,varargin)
            %Feval compute the gradient of ft
            %Input:
            %   N - total number of variables (not blocks)
            %Output:
            %   Ft(xtm1,xt) - Ft out zobjvarargineropadded to R^{nT}
            %--------  end of function help -------------%
            
            %Separate input to xtm1 and xt if input is given as one vector
            if nargin==2
                x = varargin{1};
                midInd = obj.ntm1;
                xtm1 = x(1:midInd);
                xt = x(midInd+1:end);
            else
                xtm1 = varargin{1};
                xt = varargin{2};
            end
            
            Ft =  obj.ltObj.grad(xtm1,xt) ;
        end
        
        function ft = feval (obj,varargin)
            %METHOD1
            %Input:
            %   nT - total number of variables (not blocks)
            %
            %Output:
            %   f_t(xtm1,xt) - f out zeropadded to R^{nT}
            %--------  end of function help -------------%
            if nargin==2
                x = varargin{1};
                midInd = obj.ntm1;
                xtm1 = x(1:midInd);
                xt = x(midInd+1:end);
            else
                xtm1 = varargin{1};
                xt = varargin{2};
            end
            ft =  obj.ltObj.ft(xtm1,xt);
        end
    end
end


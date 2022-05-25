classdef CellList < handle
    % forms a list array of cell objects
    
    properties (Access = private)
        
        p_eol; %point to next available link
        
    end
    
    properties
        list;
        counter;
    end
    
    methods
        function obj = CellList(K,L)
            if nargin ==0
                K = [100,1]; %default size
            end
            
            obj.list =cell(K);
            obj.p_eol = 1;
            obj.counter = 0;
            
            if nargin==2
                %                 l= numel(L);
                %                 obj.list(obj.p_eol:obj.p_eol+l)=L;
                %                 obj.p_eol = obj.p_eol + l;
                %                 obj.counter = obj.counter + l;
                obj.add(L);
            end
        end
        
        function add(obj,L)
            l= numel(L);
            obj.list(obj.p_eol:obj.p_eol+l-1)=L;
            obj.p_eol = obj.p_eol + l;
            obj.counter = obj.counter + l;
        end
        
    end
    
end


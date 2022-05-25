classdef nhpp_bspline < handle
    % untitled8 Add summary here
    %
    % NOTE: When renaming the class name untitled8, the file name
    % and constructor name must be updated to use the class name.
    %
    % This template includes most, but not all, possible properties,
    % attributes, and methods that you can implement for a System object.
    
    % Public, tunable properties
    properties
        L %spord
        support
    end
    
    % Public, non-tunable properties
    properties(Nontunable)
        
    end
    
    properties(DiscreteState)
        
    end
    
    % Pre-computed constants
    properties(Access = private)
                phi_n
    end
    
    methods
        % Constructor
        function obj = nhpp_bspline(varargin)
            % Support name-value pair arguments when constructing object
            % Perform one-time calculations, such as computing constants
            spord = varargin{1};
            bspline_L = eval(['@bspline',num2str(spord)]);  %pick bspline* where *={1,2,3,4}
            obj.L = spord;
            obj.phi_n =@(t,n) bspline_L(t-n);
            obj.support = @(n) [n, n] + [-1 , + 1 ]*(spord+1)/2; %support for one function
        end
    end
    
    methods(Access = public)
        function y = step(obj,t,n)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
            if numel(n)==1
                y = obj.stepImp(t,n);
            else
                y = obj.stepVectorImp(t,n);
            end           
        end
    end
    
    methods(Access = private)
        %% Common functions
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
        end
                
        function [outStep] = stepImp(obj,t,n)
            % filter phi for several n
                outStep= obj.phi_n(t,n);
                outStep=outStep(:).';
        end
        
        function [outMat] = stepVectorImp(obj,t,nvec)
            % filter phi for several n
            outMat = zeros(numel(nvec),numel(t));
            for i=1:numel(nvec)
                ni = nvec(i);
                outMat(i,:) = obj.phi_n(t,ni);
            end
        end
        
        
        
        function resetImpl(obj)
            % Initialize / reset discrete-state properties
        end
        
        %% Backup/restore functions
        function s = saveObjectImpl(obj)
            % Set properties in structure s to values in object obj
            
            % Set public properties and states
            s = saveObjectImpl@matlab.System(obj);
            
            % Set private and protected properties
            %s.myproperty = obj.myproperty;
        end
        
        function loadObjectImpl(obj,s,wasLocked)
            % Set properties in object obj to values in structure s
            
            % Set private and protected properties
            % obj.myproperty = s.myproperty;
            
            % Set public properties and states
            loadObjectImpl@matlab.System(obj,s,wasLocked);
        end
        
        %% Advanced functions
        function validateInputsImpl(obj,u)
            % Validate inputs to the step method at initialization
        end
        
        function validatePropertiesImpl(obj)
            % Validate related or interdependent property values
        end
        
        function ds = getDiscreteStateImpl(obj)
            % Return structure of properties with DiscreteState attribute
            ds = struct([]);
        end
        
        function processTunedPropertiesImpl(obj)
            % Perform calculations if tunable properties change while
            % system is running
        end
        
        function flag = isInputSizeLockedImpl(obj,index)
            % Return true if input size is not allowed to change while
            % system is running
            flag = true;
        end
        
        function flag = isInactivePropertyImpl(obj,prop)
            % Return false if property is visible based on object
            % configuration, for the command line and System block dialog
            flag = false;
        end
    end
    
    
    
end

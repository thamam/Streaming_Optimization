classdef MEMBUFFERCLASS<handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        maxsz
        datalist
        memCnt
        itmsz
        type
        lptr %pointer to the first(oldest) function in the active functions array
        lptrkey %holds key to which f_t is stored at datalist{lptr}
        archieve % if dbgMod = true => store all variables that were ever inserted into the buffer
        dbgMod
    end
    
    methods
        function obj = MEMBUFFERCLASS(maxsz, itmsz,itmtype,dbgMod)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            if exist('dbgMod','var')&& (~isempty(dbgMod))
                obj.dbgMod = dbgMod;
            else
                obj.dbgMod = false;
            end
            obj.maxsz = maxsz; % in units of functions 1,...,T
            obj.itmsz=itmsz;
            obj.type = itmtype;
            if (strcmp(itmtype,'double'))
                if (isinf(maxsz))
                    obj.datalist=[];
                else
                    obj.datalist = zeros(maxsz*itmsz,1);
                end
                
                obj.archieve=[];
            else
                if (isinf(maxsz))
                    obj.datalist = {};
                else
                    obj.datalist = cell(1,maxsz*itmsz);
                end
                obj.archieve={};
            end
            obj.memCnt = 0;
            obj.lptr=-1; %initialize lptrt as invalid value
        end
        
        function dout  = add(obj,din)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            assert(isa(din,obj.type));
            obj.memCnt = obj.memCnt+1;%treat each member as"one" unit
            if (strcmp(obj.type,'double'))
                if obj.dbgMod
                    obj.archieve = [obj.archieve;din];
                end
                dout = obj.datalist(1:obj.itmsz);
                tmpdatalist = circshift(obj.datalist,-obj.itmsz);
                nwWrtInd = (obj.maxsz-1)*obj.itmsz+1; %index where to start writing of new data block
                tmpdatalist(nwWrtInd:end) = din;
            else%we access a cell not mat array
                dout = obj.datalist{1:obj.itmsz};
                tmpdatalist = circshift(obj.datalist,-obj.itmsz);
                nwWrtInd = (obj.maxsz-1)*obj.itmsz+1; %index where to start writing of new data block
                tmpdatalist{nwWrtInd:end} = din;
                if obj.dbgMod
                    obj.archieve{obj.memCnt}=din;
                end
            end
            obj.datalist = tmpdatalist;
            
            if obj.memCnt<obj.maxsz
                dout=[];
            end
            if obj.lptr==(-1)
                obj.lptr=obj.maxsz;
                obj.lptrkey=obj.memCnt;
            elseif obj.lptr>1
                obj.lptr=obj.lptr-1;
            else
                %do nothing obj.lptr==1;
                obj.lptrkey=obj.maxsz-obj.memCnt+1;
            end
        end
        
        function dataout = readout(obj,ind)
            if nargin<2
                if obj.memCnt>=obj.maxsz
                    dataout = obj.datalist;
                else % Ignore empty values in readout
                    dataout = obj.datalist((obj.maxsz-obj.memCnt)*obj.itmsz+1:end);
                end
            else
                dataout = obj.datalist((1:obj.itmsz)+ind*obj.itmsz);
            end
            
            
        end
    end
end


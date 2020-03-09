classdef AtlasDrive32_1 < hardware.electrode.probe
    
    properties
    end
    
    methods
        function p = AtlasDrive32_1(varargin)
            p@hardware.electrode.probe(varargin{:}); % base class constructor
            
            p.name = 'AtlasDrive32_1';
            p.manufacturer = 'Atlas';
            p.design       = 'E32-50-S1-L6';
            p.num          = '2015-149';
            p.xcoords      = zeros(32,1);
            p.ycoords      = zeros(32,1);
            p.zcoords      = (1:32)' * 50;
            % --- specs from Atlas
%             p.channelMap   = [3 33 7 37 18 28 4 34 17 27 9 39 5 35 8 38
%             10 40 1 31 6 36 19 29 2 32 16 26 20 30 15 25];

            % --- Measured channel map
            p.channelMap = [5 21 7 23 3 19 12 28 14 30 8 24 6 22 10 26 9 25 4 20 11 27 16 32 13 29 2 18 1 17 15 31];
            p.connector    = 'MOLC';
            p.material     = 'IrOx';
            p.impedence    = 0.24;
            
            
        end
        
    end
end
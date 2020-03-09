classdef AtlasZifOmnDrive_1 < hardware.electrode.probe
    
    properties
    end
    
    methods
        function p = AtlasZifOmnDrive_1(varargin)
            p@hardware.electrode.probe(varargin{:}); % base class constructor
            
            p.name = 'AtlasZifOmnDrive_1';
            p.manufacturer = 'Atlas';
            p.design       = 'E32-50-S1-L6';
            p.num          = '2017-177';
            p.xcoords      = zeros(32,1);
            p.ycoords      = zeros(32,1);
            p.zcoords      = (1:32)' * 50;
            p.channelMap   = [18 38 23 3 17 37 24 4 16 36 25 5 15 35 26 6 14 34 27 7 13 33 28 8 12 32 29 9 30 10 11 31];
            p.connector    = 'ZIF';
            p.material     = 'IrOx';
            p.impedence    = 0.4;
            
        end
        
    end
end
classdef AtlasE32R35S1L8NT < hardware.electrode.probe
    
    properties
    end
    
    methods        
        function p = AtlasE32R35S1L8NT(varargin)
            p@hardware.electrode.probe(varargin{:}); % base class constructor
            
            p.name = 'Shank2';
            p.manufacturer = 'Atlas';
            p.design       = 'E32+R-35-S1-L8 NT';
            p.num          = '2018-154';
            p.xcoords      = zeros(32,1);
            p.ycoords      = fliplr(1:32)' * 35;
            p.zcoords      = zeros(32,1);
            p.channelMap   = [20 17 2 35 21 16 3 34 22 15 4 33 23 14 5 32 24 13 6 31 25 12 7 30 26 11 8 29 9 28 27 10];
            p.connector    = 'Omnetics36';
            p.material     = 'IrOx';
            p.impedence    = 0.5;
            
        end
        
    end
end
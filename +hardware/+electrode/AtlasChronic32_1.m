classdef AtlasChronic32_1 < hardware.electrode.probe
    
    properties
    end
    
    methods        
        function p = AtlasChronic32_1(varargin)
            p@hardware.electrode.probe(varargin{:}); % base class constructor
            
            p.name = 'AtlasChronic32_1';
            p.manufacturer = 'Atlas';
            p.design       = 'E32-50-S1-L6';
            p.num          = '2017-082';
            p.xcoords      = zeros(32,1);
            p.ycoords      = zeros(32,1);
            p.zcoords      = (1:32)' * 50;
            p.channelMap   = [20 35 21 34 22 33 23 32 24 31 25 30 26 29 27 28 2 17 3 16 4 15 5 14 9 10 8 11 6 13 7 12];
            p.connector    = 'Omnetics36';
            p.material     = 'Pt';
            p.impedence    = 0.7;
            
        end
        
    end
end
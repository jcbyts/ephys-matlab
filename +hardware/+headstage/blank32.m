classdef blank32 < hardware.headstage.headstage
    
    properties
    end
    
    methods
        function p = blank32(varargin)
            p@hardware.headstage.headstage(varargin{:}); % base class constructor
            
            p.name = 'intan_RHD2132';
            p.manufacturer = 'intan';
            p.model        = 'RHD2132';
            p.filter       = [1 7500];
            p.samplingRate = 30000;
            p.connector    = 'Omnetics36';
            p.gains        = nan;
            p.channelMap   = 1:32;
        end
        
        function chanMap = mapChannels(p, probe)
            
            assert(isa(probe, 'hardware.electrode.probe'), 'argument must be a probe')
            
            
            chanMap = p.channelMap(probe.channelMap);
            
        end
    end
end
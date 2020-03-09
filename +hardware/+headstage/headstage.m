classdef headstage < handle
    %HEADSTAGE class describes headstages
    
    properties
        name@char
        manufacturer@char
        model@char
        filter@double
        samplingRate@double
        connector@char
        gains@double
        channelMap@double
    end
    
    methods
        function h = headstage(varargin)
            
        end
        
        function chanMap = mapChannels(p, probe)
            
            assert(isa(probe,'hardware.electrode.probe'), 'mapChannels requires a probe as input')
            chanMap = p.channelMap(probe.channelMap);
            
        end
    end
    

end


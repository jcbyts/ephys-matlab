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
        
    end
    
end


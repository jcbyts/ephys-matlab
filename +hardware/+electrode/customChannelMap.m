classdef customChannelMap < hardware.electrode.probe
    % specify a custom channel map
    properties
    end
    
    methods        
        function p = customChannelMap(chanMap, varargin)
            p@hardware.electrode.probe(varargin); % base class constructor
            
            Nchan = numel(chanMap);
            
            p.name = 'customChannelMap';
            p.manufacturer = '??';
            p.design       = '??';
            p.num          = '??';
            p.xcoords      = zeros(Nchan,1);
            p.ycoords      = zeros(Nchan,1);
            p.zcoords      = zeros(Nchan,1);
            p.channelMap   = chanMap;
            p.connector    = '??';
            p.material     = '??';
            p.impedence    = nan;
            p.headstages   = {hardware.headstage.blank32};
            
        end
        
    end
end
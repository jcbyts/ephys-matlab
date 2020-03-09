classdef SingleElectrodes < hardware.electrode.probe
    
    properties
    end
    
    methods        
        function p = SingleElectrodes(chanMap, varargin)
            p@hardware.electrode.probe(varargin); % base class constructor
            
            Nchan = numel(chanMap);
            
            p.name = 'SingleElectrodes';
            p.manufacturer = '??';
            p.design       = '??';
            p.num          = '??';
            p.xcoords      = zeros(Nchan,1);
            p.ycoords      = zeros(Nchan,1);
            p.zcoords      = zeros(Nchan,1);
            p.channelMap   = chanMap;
            p.connector    = 'Omnetics36';
            p.material     = '??';
            p.impedence    = nan;
            p.headstages   = hardware.headstage.blank32;
            
        end
        
    end
end
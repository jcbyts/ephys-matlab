classdef intan_RHD2132 < hardware.headstage.headstage
    
    properties
    end
    
    methods
        function p = intan_RHD2132(varargin)
            p@hardware.headstage.headstage(varargin{:}); % base class constructor
            
            p.name = 'intan_RHD2132';
            p.manufacturer = 'intan';
            p.model        = 'RHD2132';
            p.filter       = [1 7500];
            p.samplingRate = 30000;
            p.connector    = 'Omnetics36';
            p.gains        = nan;
            p.channelMap   = [0 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 0 0 25 26 27 28 29 30 31 32 1 2 3 4 5 6 7 8 0];
        end
        
         function chanMap = mapChannels(p, probe)
            chanMap = [];
            
            assert(isa(probe, 'hardware.electrode.probe'), 'argument must be a probe')
            
            switch probe.connector
                
                case 'ZIF'
                    
                    Zif2Omn = [0 0 35 34 33 32 31 30 29 28 27 26 25 24 23 22 21 20 00 00 00 00 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17];
                    chanMap = p.channelMap(Zif2Omn(probe.channelMap));
                    
                case 'Omnetics36'
                    
                    chanMap = p.channelMap(probe.channelMap);
                    
                case 'MOLC'
                    warning('Assuming this is a measured channel map. No transformation applied')
                    chanMap = probe.channelMap;
                    
            end
            
            
        end
    end
end
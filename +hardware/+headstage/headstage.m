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
        
        %         function chanMap = channelMap(h, probe)
        %             chanMap = [];
        %
        %             assert(isa(probe, @hardware.electrode.probe), 'argument must be a probe')
        %
        %             switch probe.connector
        %
        %                 case 'ZIF'
        %
        %                 case 'Omnetics36'
        %
        %                 case 'MOLC'
        %
        %             end
        %
        %             switch h.connector
        %
        %             end
        %
        %
        %         end
        %
        %         handles.ops.chanMap = [5 21 7 23 3 19 12 28 14 30 8 24 6 22 10 26 9 25 4 20 11 27 16 32 13 29 2 18 1 17 15 31];
        %     case 'Atlas32ChZIF'
        %         handles.ops.chanMap = [8 24 9 25 7 23 10 26 6 22 11 27 5 21 12 28 4 20 13 29 3 19 14 30 2 18 15 31 16 32 1 17];
        %
        %         Zif = [18 38 23 3 17 37 24 4 16 36 25 5 15 35 26 6 14 34 27 7 13 33 28 8 12 32 29 9 30 10 11 31];
        %         Zif2Omn = [0 0 35 34 33 32 31 30 29 28 27 26 25 24 23 22 21 20 00 00 00 00 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17];
        %         intanHeadstage = [0 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 0 0 25 26 27 28 29 30 31 32 1 2 3 4 5 6 7 8 0];
        %
        %         handles.ops.chanMap = intanHeadstage(Zif2Omn(Zif));
        %
        %
        %     case 'Atlas32ChOmnChronic'
        %
        %         rawChanMap = [20 35 21 34 22 33 23 32 24 31 25 30 26 29 27 28 2 17 3 16 4 15 5 14 9 10 8 11 6 13 7 12];
        %         intanHeadstage = [0 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 0 0 25 26 27 28 29 30 31 32 1 2 3 4 5 6 7 8 0];
        %
        %         handles.ops.chanMap = intanHeadstage(rawChanMap);
        %
        %     otherwise
        %         handles.ops.chanMap = 1:handles.ops.Nchan;
        %     end
    end
    
end


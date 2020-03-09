classdef probe < handle
    % PROBE  Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name@char
        manufacturer@char
        design@char
        num@char
        xcoords@double
        ycoords@double
        zcoords@double
        channelMap@double
        connector@char
        material@char
        impedence@double
        headstages@cell
    end
    
    methods
        
        function p = probe(varargin)
           
        end
    end
    
end


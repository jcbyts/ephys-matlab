classdef eyepos
    % Class for common operations on eye tracking data.
    %
    % To see the public properties of this class, type
    %
    %   properties(marmodata.eye)
    %
    % To see a list of methods, type
    %
    %   methods(marmodata.eye)
    %
    % The class constructor may be called with a range of arguments:
    %
    %   d = marmodata.eye;                       % an empty @eye object.
    %   d = marmadata.eye(STRUCT);               % an @eye object initialized with the contents of STRUCT
    %   d = marmodata.eye(TIME,X,Y,PWDTH,PHGHT); % an @eye object initialised with gaze position (X,Y) and pupil dimensions (PWDTH,PHGHT)
    
    % 2016-11-14 - Shaun L. Cloherty <s.cloherty@ieee.org>
    
    properties (SetAccess = private, GetAccess = public)
        tsample@double; % raw sample times (no offset/aligment)
    end
    
    properties
        %     t@double; % time
        
        x@double; % gaze position (horiz.)
        y@double; % gaze position (vert.)
        
        pwdth@double; % pupil width
        phght@double; % pupil height
        
        toffset@double = 0; % time of first sample
    end
    
    % dependent properties, calculated on the fly...
    properties (Dependent, SetAccess = private, GetAccess = public)
        t; % time
        
        fs; % sampling freq. (samples/s; Hz), fs = 1/dt
        dt; % sampling interval (s/sample), dt = 1/fs
        parea; % pupil area
    end
    
    methods
        function value = get.t(d)
            value = d.tsample + d.toffset;
        end
    end
    
    methods % get/set dependent properties
        % dependent property get methods
        function value = get.fs(d)
            % sampling frequency (samples/s)
            value = 1./d.dt;
        end
        
        function value = get.dt(d)
            % sampling interval (s/sample)
            value = median(diff(d.t));
        end
        
        function value = get.parea(d)
            % pupil area (a.u.)
            value = pi*d.pwdth.*d.phght; % area of an ellipse
        end
    end
    
    methods (Access = public)
        
        function d = eyepos(varargin) % constructor
            if nargin == 0
                return
            end
            
            switch nargin
                case 1
                    d_ = varargin{1};
                    if ~isstruct(d_)
                        error('Was expecting a struct but got %s instead!',class(d_));
                    end
                    for fname = {'dt','x','y','pwdth','phght','toffset'}
                        d.(fname{1}) = d_.(fname{1});
                    end
                case 5
                    d.tsample = varargin{1}; % time (seconds)
                    
                    d.x = varargin{2}; % gaze position
                    d.y = varargin{3};
                    
                    d.pwdth = varargin{4}; % pupil wdth
                    d.phght = varargin{5}; % pupil hght
                otherwise
                    error('Invailid arguments in call to marmodata.eyedata().');
            end
        end
        
        function value = t2s(d,t)
            % time to sample/index
            [~,value] = min(abs(d.t - t));
        end
        
        function value = s2t(d,s)
            % sample/index to time
            value = d.t(s);
        end
        
    end % public methods
    
    methods (Static)
        
        function obj = loadobj(s)
            if isstruct(s)
                % we end up here when loading 'old' eye objects
                % that don't contain tsample
                assert(~isfield(s,'tsample'),'Error loading @eye object...');
                
                % if s contains toffset then s.t ws saved with this offset
                % applied... therefore, we need to subtract toffset
                toffset = 0;
                if isfield(s,'toffset')
                    toffset = s.toffset;
                end
                t = s.t - toffset;
                
                obj = marmodata.eye(t,s.x,s.y,s.pwdth,s.phght);
                obj.toffset = toffset;
            else
                obj = s;
            end
        end
        
    end % static methods
    
end % classdef
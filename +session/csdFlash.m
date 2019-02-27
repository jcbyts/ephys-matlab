classdef csdFlash < handle
    %CSDFLASH import module for csd flash protocol
    %   csdFlash imports the csd flash stimulus module into a common format
    %   It can also analyze the csd and find sinks and sources
    
    properties
        numTrials
        trial
        display
        method@char='spline' % csd method
        sessionTrialIdx=[]
    end
    
    methods
        
        function obj = csdFlash(PDS)
            
            % --- find CSD flash trials
            stim = 'csdFlash';
            
            [hasStim, numTrialsPerPDS] = io.findPDScontainingStimModule(PDS, stim);
            trialOffset = [0; cumsum(numTrialsPerPDS)];
            
            if ~any(hasStim)
                return
            end
            
            obj.display = PDS{find(hasStim, 1, 'first')}.initialParametersMerged.display;
            
            
            for i = find(hasStim(:)')
                
                [trial_, ~, idx] = obj.importPDS(PDS{i});
                
                obj.sessionTrialIdx = [obj.sessionTrialIdx trialOffset(i) + idx(:)'];
                
                if isempty(trial_)
                    continue
                end
                
                obj.trial = [obj.trial; trial_(:)];
                
            end
            
            
            obj.numTrials = numel(obj.trial);
            
        end
        
        function stats = computeCsd(obj, ops, varargin)
            % takes in ops
%              stats = computeCsd(obj, ops, varargin)
            ip = inputParser();
            ip.addOptional('plot', false)
            ip.parse(varargin{:})
            
            [lfp, ~, lfpInfo] = io.getLFP(ops);
            
            flashTimes = [obj.trial.onset];
            
            chmap = load(ops.chanMap);
            
            stats = obj.csdBasic(lfp, flashTimes, lfpInfo, 'method', obj.method, 'channelDepths', chmap.ycoords, 'plot', ip.Results.plot);
            
        end
    end
    
    methods (Static)
        function [trial,display,trialIdx] = importPDS(PDS)
            
            pdsDate = PDS.initialParametersMerged.session.initTime;
            
            if pdsDate > datenum(2018,02,01)
                [trial,display,trialIdx] = session.csdFlash.importPDS_v2(PDS);
            else
                [trial,display,trialIdx] = session.csdFlash.importPDS_v1(PDS);
            end
                
        end
        
        function [csdTrial,display,stimTrials] = importPDS_v2(PDS)
            csdTrial = struct();
            
            % --- find CSD flash trials
            stim = 'csdFlash';
            
            trialIx = cellfun(@(x) isfield(x, stim), PDS.data);
            display = [];
            stimTrials = find(trialIx);
            
            for j = 1:numel(stimTrials)
                thisTrial = stimTrials(j);
                
                kTrial = j;
                
                csdTrial(kTrial).frameTimes = PDS.PTB2OE(PDS.data{thisTrial}.timing.flipTimes(1,1:end-1));
                csdTrial(kTrial).start      = csdTrial(kTrial).frameTimes(1);
                csdTrial(kTrial).stop       = csdTrial(kTrial).frameTimes(end);
                csdTrial(kTrial).duration   = PDS.PTB2OE(PDS.data{thisTrial}.timing.flipTimes(1,end-1)) - csdTrial(kTrial).start;
                
                if isempty(PDS.data{thisTrial}.(stim).hFlash.log)
                    csdTrial(kTrial).onset = [];
                    csdTrial(kTrial).offset = [];
                    continue
                end
                
                onIx = PDS.data{thisTrial}.(stim).hFlash.log(1,:)==1;
                offIx = PDS.data{thisTrial}.(stim).hFlash.log(1,:)==0;
                offIx(1) = false;
%                 csdTrial(kTrial).on         = PDS.data{thisTrial}.(stim).on;
                
                % TODO: Align to frame onset
                csdTrial(kTrial).onset      = PDS.PTB2OE(PDS.data{thisTrial}.(stim).hFlash.log(2,onIx));
                csdTrial(kTrial).offset     = PDS.PTB2OE(PDS.data{thisTrial}.(stim).hFlash.log(2,offIx));
            end
            
            
        end
           
        function [csdTrial,display,stimTrials] = importPDS_v1(PDS)
            csdTrial = struct();
            
            % --- find CSD flash trials
            stim = 'csdFlash';
            
            trialIx = cellfun(@(x) isfield(x, stim), PDS.data);
            display = [];
            stimTrials = find(trialIx);
            
            for j = 1:numel(stimTrials)
                thisTrial = stimTrials(j);
                
                kTrial = j;
                
                csdTrial(kTrial).frameTimes = PDS.PTB2OE(PDS.data{thisTrial}.timing.flipTimes(1,1:end-1));
                csdTrial(kTrial).start      = csdTrial(kTrial).frameTimes(1);
                csdTrial(kTrial).stop       = csdTrial(kTrial).frameTimes(end);
                csdTrial(kTrial).duration   = PDS.PTB2OE(PDS.data{thisTrial}.timing.flipTimes(1,end-1)) - csdTrial(kTrial).start;
                
                csdTrial(kTrial).on         = PDS.data{thisTrial}.(stim).on;
                csdTrial(kTrial).onset      = csdTrial(kTrial).frameTimes(diff(csdTrial(kTrial).on)==1);
            end
            
            
        end
        
        function stats = csdBasic(lfp, eventTimes, lfpInfo, varargin)
            % CSDBASIC computes the current source density
            % Inputs:
            %   LFP     [nTime x nChannel] - raw voltage traces
            %   events       [nEvents x 1] - timestamps of the events
            %   lfpInfo           [struct] - lfp info struct from io.getLfp
            %
            % optional arguments (as argument pairs):
            %   'channelDepths'  [nChannels x 1] - array of channel depths
            %   'window'         [1 x 2]         - start and stop of analysis window
            %                                      (aligned to event time)
            %   'plot'           logical         -  plot if (default: true)
            %   'method'         string          - csd method (default: 'spline')
            %
            % valid csd methods:
            %       'standard' - second spatial derivative
            %       'step'     - stepwise inverse method (not really sure)
            %       'spline'   - interpolated inverse CSD
            %
            % 2017 jly wrote it
            
            ip = inputParser();
            ip.addParameter('window', [-100 200])
            ip.addParameter('channelDepths', [])
            ip.addParameter('plot', true)
            ip.addParameter('method', 'spline')
            ip.parse(varargin{:});
            
            eventTimes = eventTimes(:);
            if ~all(mod(lfpInfo.fragments(:),1)==0)
                warning('Rounding LFP fragments. This is a hack')
                lfpInfo.fragments = round(lfpInfo.fragments); % This is a hack
            end
            % conver times to samples
            ev = io.convertTimeToSamples(eventTimes, lfpInfo.sampleRate, lfpInfo.timestamps(:), lfpInfo.fragments(:));
            
            % event-triggered LFP
            [sta,~, time] = pdsa.eventTriggeredAverage(lfp, ev(:), ip.Results.window);
            
            if isempty(ip.Results.channelDepths)
                ch0 = (-32:-1)*50;
            else
                ch0 = ip.Results.channelDepths(:)'; % row vector
            end
            
            switch ip.Results.method
                case 'spline'
                    % compute the CSD using the spline inverse method
                    CSD = csd.splineCSD(sta', 'el_pos', ch0);
                case 'standard'
                    CSD = csd.standardCSD(sta', 'el_pos', ch0);
                case 'step'
                    CSD = csd.stepCSD(sta', 'el_pos', ch0);
                otherwise
                    error('valid methods are {spline, standard, step}')
            end
            
            % find the sink and reversal point
            ix = time > 0 & time < 100; % look over time window after flash
            
            % sink should be the minimum value
            [~,id] = min(reshape(CSD(:,ix), [], 1));
            % convert to indices
            [depthIndex,timeIndex] = ind2sub(size(CSD(:,ix)), id);
            % upsample channels to index into them
            chUp   = linspace(1, numel(ch0), size(CSD,1));
            depthUp= linspace(ch0(1), ch0(end), size(CSD,1));
            
            % find reversal point
            CSD_ = CSD(:,ix);
            reversalPoints = findZeroCrossings(CSD_(:,timeIndex));
            
            % output structure
            stats.STA   = sta';
            stats.CSD   = CSD;
            stats.time  = time;
            stats.depth = depthUp;
            stats.chDepths = ch0;
            stats.chUp  = chUp;
            stats.sinkDepth = depthUp(depthIndex);
            stats.sinkChannel = chUp(depthIndex);
            stats.reversalPointDepth = depthUp(reversalPoints);
            
            if ip.Results.plot % afterall, it is a plot function
                imagesc(time, ch0, CSD-mean(CSD(:))); axis xy
                colormap jet
                hold on
                plot(time, bsxfun(@plus, sta, ch0), 'Color', repmat(.1, 1, 3))
                xlim(ip.Results.window)
                plot(time([1 end]), stats.sinkDepth*[1 1], 'w--', 'Linewidth', 2)
                plot(time([1 end]), [1; 1]*stats.reversalPointDepth, 'r--', 'Linewidth', 2)
                tmp = abs(stats.reversalPointDepth - stats.sinkDepth);
                tmp = tmp + stats.sinkDepth;
                plot(time([1 end]), [1; 1]*tmp, 'r--', 'Linewidth', 2)
                axis ij
            end
        end

        function i = findZeroCrossings(data, mode)
            %FINDZEROCROSSINGS Find zero crossing points.
            %   I = FINDZEROCROSSINGS(DATA,MODE) returns the indicies into the supplied
            %   DATA vector, corresponding to the zero crossings.
            %
            %   MODE specifies the type of crossing required:
            %     MODE < 0 - results in indicies for the -ve going zero crossings,
            %     MODE = 0 - results in indicies for ALL zero crossings (default), and
            %     MODE > 0 - results in indicies for the +ve going zero crossings.
            
            % $Id: findZeroCrossings.m,v 1.1 2008-07-21 23:31:50 shaunc Exp $
            
            if nargin < 2
                mode = 0;
            end
            
            [i,~,p] = find(data); % ignore zeros in the data vector
            
            switch sign(mode)
                case -1
                    % find -ve going crossings
                    ii = find(diff(sign(p))==-2);
                case 0
                    % find all zero crossings
                    ii = find(abs(diff(sign(p)))==2);
                case 1
                    % find +ve going crossings
                    ii = find(diff(sign(p))==2);
            end;
            
            i = round((i(ii)+i(ii+1))/2);
        end
    end
end

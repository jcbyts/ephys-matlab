function [data, timestamps, elInfo] = getVpx(sess, PDS, overwrite)
% LOAD Load data from a ViewPoint data file.
% The first time it is called, getVpx converts the vpx file to a binary/mat
% combo, where the binary file contains raw data and the mat file contains
% meta-data. 
% Inputs:
%   session@struct    - session info struct (or struct-array)
%   PDS@cell          - cell array of PDS structures
%   overwrite@logical - flag to reimport
% Outputs:
%   data@double       - n x 3 [x y pupil] in degrees
%   timestamps@double - n x 1 (time in OE time)
%   elInfo@struct     - meta data
% Example call:
% [data, timestamps, elInfo] = io.getVpx(ops, PDS)

% 2017.08.14     jly     wrote it


if nargin < 3
    overwrite = false;
end

flipX = false;
flipY = false;

if isa(sess, 'table') % it's meta data -- load the struct
    if ~isnan(sess.FlipXEye) && sess.FlipXEye
        flipX = true;
    end
    
    if ~isnan(sess.FlipYEye) && sess.FlipYEye
        flipY = true;
    end
    
    sess = io.loadSession(sess);
end

% initialize info struct
elInfo.timestamps = [];
elInfo.fragments  = [];
elInfo.sampleRate = [];
elInfo.dateNum    = [];
elInfo.fields     = {'EyeX', 'EyeY', 'PupilArea'};
elInfo.bitDeg     = 1e-3;

feye    = fullfile(sess.path, '_behavior', 'eyepos.dat');
felinfo = fullfile(sess.path, '_behavior', 'eye_info.mat');

% if the import has already been run, just load the existing data
if exist(feye, 'file') && ~overwrite
    elInfo = load(felinfo);
    
    fid = fopen(feye, 'r');
    fseek(fid, 0, 'eof');
    filesize = ftell(fid);
    fseek(fid, 0, 'bof');
    
    nFields = numel(elInfo.fields);
    nTotSamps = filesize/nFields/2;
    
    data = double(fread(fid, [nFields nTotSamps], '*uint16'));
    
    [data, timestamps] = convertRawToDeg(data, elInfo, flipX, flipY);
    
    return
    
end

if nargin < 2
    PDS = io.getPds(sess);
end

fidout = fopen(feye, 'w');

nPds = numel(PDS);

for kPds = 1:nPds
    
%     try
    [dataRAW, ~, el_info] = getVpxData(sess, PDS{kPds});
    fwrite(fidout, dataRAW', '*uint16') % write to disk
    
    % update info struct
    elInfo.timestamps   = [elInfo.timestamps el_info.timestamps];
    elInfo.fragments    = [elInfo.fragments el_info.fragments];
    elInfo.sampleRate   = [elInfo.sampleRate el_info.sampleRate];
    elInfo.dateNum      = [elInfo.dateNum el_info.dateNum];
    elInfo.bitDeg       = el_info.bitDeg;
%     catch
%         warning('no edf data for PDS %d', kPds)
%     end
    
end

elInfo.sampleRate = unique(elInfo.sampleRate);
%assert(numel(elInfo.sampleRate)==1, 'Eyelink sample rate changed throughout session. This is not supported')

fclose(fidout);
save(felinfo, '-v7.3', '-struct', 'elInfo');

fid = fopen(feye, 'r');
fseek(fid, 0, 'eof');
filesize = ftell(fid);
fseek(fid, 0, 'bof');

nFields = numel(elInfo.fields);
nTotSamps = filesize/nFields/2;

data = double(fread(fid, [nFields nTotSamps], '*uint16'));
[data, timestamps] = convertRawToDeg(data, elInfo, flipX, flipY);

% -------------------------------------------------------------------------
% Helper functions
%--------------------------------------------------------------------------

% --- Convert from raw file values to degrees and time
function [data, timestamps] = convertRawToDeg(data, elInfo, flipX, flipY)
nTotSamps = size(data,2);
data(1:2,:) = data(1:2,:)*elInfo.bitDeg(1) + elInfo.bitDeg(2);
if ~all(data(3,:)==0)
    data(1,data(3,:)==0) = nan;
    data(2,data(3,:)==0) = nan;
end

timestamps = io.convertSamplesToTime(1:nTotSamps, elInfo.sampleRate(1), elInfo.timestamps(:), elInfo.fragments(:));

if flipX
    data(1,:) = -data(1,:);
end

if flipY
    data(2,:) = -data(2,:);
end

% --- Apply the calibration matrix to data read in from the vpx file
function [data, timestamps, info] = getVpxData(sess, PDS)
% INPUT
%   sess
%   PDS
%
% OUTPUT
%   data       [n x 3] - the raw data [X, Y, Pupil]
%   timestamps [n x 1] - the time of each sample
%   info      [struct] - info (sampling rate, timestamps, fragments)
%
% At the top of each Vpx file is file header information that includes
% the date and time that the data was collected, the apparatus geometry
% settings that can be used to obtain the gaze angle in degrees, whether
% smoothed or unsmoothed data was stored, etc.
%
% Each line of the data file is a unique data record. The type of record
% is indicated by an integer tag value in the first column. The tag
% determines how the record entries on that line should be interpreted.
% The record entries on the line are tab delimited.
%
% The various data records are described in Section 13.3.4 of the
% ViewPoint user guide.

% Find all times the calibration matrix changed
if ~isfield(PDS, 'initialParameters')
    PDS.initialParameters = PDS.data;
    PDS.initialParameterNames = cell(1,numel(PDS.data));
    for i = 1:numel(PDS.data)
        PDS.initialParameterNames{i} = ['Pause after ' num2str(i)];
    end
end
arringtonChanged = find(cellfun(@(x) isfield(x, 'arrington'), PDS.initialParameters));
pauseTrials    = find(cellfun(@(x) ~isempty(regexp(x, 'Pause', 'once')), PDS.initialParameterNames));
arringtonChangedAfterPause = intersect(arringtonChanged, pauseTrials);
changedParamNames = PDS.initialParameterNames(arringtonChangedAfterPause);
trialNums = cellfun(@(x) str2double(cell2mat(regexp(x, '\d+', 'match'))), changedParamNames);

calibMatChanged = cellfun(@(x) isfield(x.arrington, 'calibration_matrix'), PDS.initialParameters(arringtonChangedAfterPause));


calibMat = [];
if ~isempty(PDS.initialParametersMerged.arrington.calibration_matrix)
    calibMat = {PDS.initialParametersMerged.arrington.calibration_matrix};
end

calibMat = [calibMat cellfun(@(x) x.arrington.calibration_matrix,...
            PDS.initialParameters(arringtonChangedAfterPause), 'UniformOutput', false)];

calibMatChangeIdx = trialNums(calibMatChanged);


VpxTimeCmChanged = [0 cellfun(@(x) x.timing.arringtonStartTime(2), PDS.data(calibMatChangeIdx))]*1e3;
%VpxTimeCmChanged = [0 cellfun(@(x) x.timing.arringtonStartTime(1), PDS.data(calibMatChangeIdx))]*1e3;

eyeIdx = PDS.initialParametersMerged.arrington.eyeIdx;


[~, fname, ~] = fileparts(PDS.initialParametersMerged.session.file);

vpxFile=fullfile(sess.path, [fname '.vpx']);
hasEdf=exist(vpxFile, 'file');
if hasEdf
    v=importVpx(vpxFile);
    vpxTrialStarts=cellfun(@(x) strcmp(x.eventId, 'TRIALSTART'), v.markers);
    vpxTrialTS=double(cellfun(@(x) x.time, v.markers(vpxTrialStarts)));
end

pdsTrialStarts=cellfun(@(x) x.timing.arringtonStartTime(1), PDS.data);

vpxTrialTS(:) = vpxTrialTS(:)*1000;  % go from secs to ms

AR2OEfit =[vpxTrialTS(:) ones(numel(vpxTrialTS),1)]\PDS.PTB2OE(pdsTrialStarts(:));
AR2OE  =@(x) x*AR2OEfit(1) + AR2OEfit(2); % this converts from Eyelink sample times to OE times

vpxTime=double(v.time);

x=double(v.x);
y=double(v.y);

numCmChanges = numel(calibMat);

% build raw input
X=[x(:)'; y(:)'; ones(1,numel(y))];

eyeXyDeg = zeros(size(X,2), 2);

for iCalib = 1:numCmChanges
    
    cm=calibMat{iCalib}(:,:,eyeIdx);
    
    if iCalib == numCmChanges
        iix = vpxTime > VpxTimeCmChanged(iCalib);
    else
        iix = vpxTime > VpxTimeCmChanged(iCalib) & vpxTime < VpxTimeCmChanged(iCalib+1);
    end
    
    eyePxTmp = cm*X(:,iix); % convert raw to pixels
    
    win = PDS.initialParametersMerged.display.winRect(3:4);
    % nan times when the eye was off the screen
    offScreen = (eyePxTmp(1,:)>=win(1) | eyePxTmp(1,:)<=0) | (eyePxTmp(2,:)>=win(2) | eyePxTmp(2,:)<=0);
    
    % center
    eyePxTmp = bsxfun(@minus, eyePxTmp, PDS.initialParametersMerged.display.ctr(1:2)');
    
    % convert to degrees
    eyeDegTmp = pds.px2deg(eyePxTmp, PDS.initialParametersMerged.display.viewdist, PDS.initialParametersMerged.display.px2w);
    eyeDegTmp(:,offScreen) = nan;
    
    eyeXyDeg(iix, :) = eyeDegTmp';
end

dT = diff(vpxTime);
ss = mode(dT);

info.sampleRate = 1e3/ss;

tmp = pds.px2deg(-PDS.initialParametersMerged.display.ctr(1:2)', PDS.initialParametersMerged.display.viewdist, PDS.initialParametersMerged.display.px2w);
tmp = min(tmp);

info.bitDeg = [1/100 tmp];

% flip Y so up is positive
eyeXyDeg(:,2) = eyeXyDeg(:,2)*-1;

data = [int16((eyeXyDeg - tmp)*100) int16((v.pwdth.*v.phght)*100)];

tol = .5; % half a millisecond tolerance
breaks = find( abs(dT - ss) > tol)'+1;

info.fragments = diff([0 breaks numel(vpxTime)]);
timestamps = AR2OE(vpxTime(:)');
info.timestamps = timestamps([1 breaks]);
info.dateNum = PDS.initialParametersMerged.session.initTime;


% --- Read in raw data fromthe Vpx File
function v = importVpx(vpxFile, varargin)
% import eye data from the Vpx File

% parse input
p = inputParser;
p.addParameter('verbose',false,@islogical);
p.addParameter('loadSamples',true,@islogical);
p.addParameter('loadMarkers',true,@islogical);
p.parse(varargin{:});

args = p.Results;

assert(exist(vpxFile,'file')==2, sprintf('%s does not exist.',vpxFile))

if args.verbose, disp('Loading data from .vpx file.'); end

try
    fid = fopen(vpxFile,'r');
catch
    error('Could not open %s.', vpxFile);
end

% colums to retrieve:
%
% 'ATT' - Eye A TotalTime
% 'ADT' - Eye Delta Time - higher precision (time since last sample)
% 'ALX' - Eye A Gaze X
% 'ALY' - Eye A Gaze Y
% 'APW' - Eye A, Pupil Width
% 'APH' - Eye A, Pupil Height
cols = {'ATT', 'ADT', 'ALX','ALY','APW','APH'};
idx = [];
      
samples = {}; % the eye data (column order as in cols above)
markers = {}; % the external event markers
      
while ~feof(fid)
  txt = fgetl(fid);
        
  % extract the data record tag
  pat = '(?<tag>[\d]+)\t+(?<payload>.*)';
  tokens = regexp(txt,pat,'names');
        
  switch tokens.tag
    case '10'
      % eye data
      pat = '[-\d\.]+';
      matches = regexp(tokens.payload,pat,'match');
      
      samples{end+1,1} = cellfun(@(x) str2double(x), matches(idx));
    case '12'
      % async. string data <-- i.e., marmoview event markers
            
      % parcel out the marmoview event markers
      %
      % an event marker is a string containing an eventId optionally
      % followed by one or more name-value parameter pairs, i.e.,
      %
      %   <eventId>[[:<param1>:<value1>] ... :<paramN>:<valueN>]]
      %
      % for example:
      %
      % TRIALSTART:TRIALNO:1
      
      pat = '(?<time>[\d\.]+)\t(?<eventId>\w+):*(?<payload>.*)';
      marker = regexp(tokens.payload,pat,'names');
      
      % break out the payload into a struct
      tmp = regexp(marker.payload,':','split');
      tmp = tmp(find(cellfun(@(x) ~isempty(x),tmp)));
      
      if mod(length(tmp),2) == 1 % length(tmp) is odd!
        warning('Malformed event payload!');
        continue;
      end

      marker.time = str2double(marker.time);
      marker.payload = struct(tmp{:});
      markers{end+1,1} = marker;
    case '2'
      % async. marker data
    case '3'
      % header
    case '5'
      % column headings
    case '6'
      % column heading abbrev.
      abbrev = regexp(tokens.payload,'\t','split');
      abbrev = abbrev(find(cellfun(@(x) ~isempty(x),abbrev)));
      idx = cellfun(@(x) find(strcmp(x,abbrev)),cols);
    case '7'
      % Unknown!?
    case '14'
      % async. head tracker data
    case '16'
      % picture image file name
    case '777'
      % movie frame number
	otherwise
      error('Unknown tag %s.',tokens.tag);
  end
end

fclose(fid);

samples = cell2mat(samples);

v.time = cumsum(samples(:,2));

v.x = samples(:,3);
v.y = samples(:,4);

v.pwdth = samples(:,5);
v.phght = samples(:,6);

v.markers = markers;

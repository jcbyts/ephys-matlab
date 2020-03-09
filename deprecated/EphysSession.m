function varargout = EphysSession(varargin)
% EPHYSSESSION MATLAB code for EphysSession.fig
%      EPHYSSESSION, by itself, creates a new EPHYSSESSION or raises the existing
%      singleton*.
%
%      H = EPHYSSESSION returns the handle to a new EPHYSSESSION or the handle to
%      the existing singleton*.
%
%      EPHYSSESSION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EPHYSSESSION.M with the given input arguments.
%
%      EPHYSSESSION('Property','Value',...) creates a new EPHYSSESSION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EphysSession_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EphysSession_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EphysSession

% Last Modified by GUIDE v2.5 24-Jul-2017 08:51:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EphysSession_OpeningFcn, ...
                   'gui_OutputFcn',  @EphysSession_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before EphysSession is made visible.
function EphysSession_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EphysSession (see VARARGIN)

% Choose default command line output for EphysSession
handles.output = hObject;
if nargin > 3
    handles.sessionDir = varargin{1};
end

% --- Turn off panels until data is loaded
handles.Panel_ArtRem.Visible    = 'off';
handles.Panel_ChanMap.Visible   = 'off';
handles.Panel_Plot.Visible      = 'off';
handles.Panel_Reference.Visible = 'off';

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes EphysSession wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = EphysSession_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Read from buffer
function plot_buffer = getBuffer(handles)
    
% draw
fseek(handles.plot_fid, handles.fid_position, -1);

plot_buffer = double(fread(handles.plot_fid, handles.plot_bufferSize, '*int16'));

[plot_buffer, handles.ops] = io.preprocessVops(plot_buffer, handles.ops);


% --- Update the traces
function handles = drawTraces(handles)

handles.plot_buffer     = getBuffer(handles);

% for i = 1:handles.ops.Nchan
%     plot_buffer(i,:) = plot_buffer(i,:)*handles.plot_yscale + handles.plot_yoffset*i;
% end

for i = 1:handles.ops.Nchan
    handles.plot_handle(i).XData = 1:numel(handles.plot_buffer(i,:));
    handles.plot_handle(i).YData = handles.plot_buffer(i,:)*handles.plot_yscale + handles.plot_yoffset*i;
end

handles.axes1.XLim = [1 numel(handles.plot_buffer(1,:))];

handles.scaleBarX = .9*diff(handles.axes1.XLim) + handles.axes1.XLim(1)*[1 1];
handles.scaleBarY = .05*diff(handles.axes1.YLim) + [0 200*handles.plot_yscale];
handles.scaleBarHandle.XData = handles.scaleBarX;
handles.scaleBarHandle.YData = handles.scaleBarY;
handles.scaleBarText.Position(1) = handles.scaleBarX(1);
handles.scaleBarText.Position(2) = mean(handles.scaleBarY);

disp(handles.fid_position)

% handles.Slider_Vertical.Max = max(handles.plot_buffer(:)) + handles.plot_yoffset;
% handles.Slider_Vertical.Min = -100;

% handles.Slider_Horizontal.SliderStep = [handles.plot_winWidth/handles.fid_length/100 handles.plot_winWidth/handles.fid_length/100];


% --- Executes on slider movement.
function Slider_Horizontal_Callback(hObject, eventdata, handles)
% hObject    handle to Slider_Horizontal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.fid_position = ceil(handles.Slider_Horizontal.Value)*handles.ops.Nchan;

handles = drawTraces(handles);

guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function Slider_Horizontal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Slider_Horizontal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function Slider_Vertical_Callback(hObject, eventdata, handles)
% hObject    handle to Slider_Vertical (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = handles.Slider_Vertical.Value;
delta = val - handles.Ypos;
handles.axes1.YLim = handles.axes1.YLim + sign(delta)*handles.plot_yoffset/2;
handles.Ypos = val;
% handles.axes1.YLim =;
% handles.plot_yoffset
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function Slider_Vertical_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Slider_Vertical (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function Input_SessionDirectory_Callback(hObject, eventdata, handles)
% hObject    handle to text_SessionDirectory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_SessionDirectory as text
%        str2double(get(hObject,'String')) returns contents of text_SessionDirectory as a double


% --- Executes during object creation, after setting all properties.
function text_SessionDirectory_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_SessionDirectory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Button_BrowseForSessionDirectory.
function Button_BrowseForSessionDirectory_Callback(hObject, eventdata, handles)
% hObject    handle to Button_BrowseForSessionDirectory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.sessionDir = uigetdir();
[~, handles.sessionName] = fileparts(handles.sessionDir);
handles.text_SessionDirectory.String = handles.sessionDir;

handles.DirectoryTitle.String = 'Loading/Converting to binary';
guidata(hObject, handles)

if exist(fullfile(handles.sessionDir, 'ops.mat'), 'file')
    handles.ops = load(fullfile(handles.sessionDir, 'ops.mat'));
else
    handles.ops        = buildConfigFile(handles.sessionDir);
end
handles.DirectoryTitle.String = handles.sessionName;

% read in binary data
handles.plot_fid        = fopen(handles.ops.fbinary);
handles.fid_position    = ftell(handles.plot_fid);
fseek(handles.plot_fid, 0,'eof');
handles.fid_length      = ftell(handles.plot_fid);
fseek(handles.plot_fid, handles.fid_position, 'bof');
handles.plot_winWidth   = 10e3;
handles.plot_bufferSize = [handles.ops.Nchan handles.plot_winWidth];
handles.plot_yscale     = 1;
handles.plot_yoffset    = 100;
fseek(handles.plot_fid, handles.fid_position, -1);

% --- default artifact removal
handles.ops.artifactThresh = inf;
handles.ops.artifactNchans = handles.ops.Nchan;

% --- Default Channel Map
fl = dir('+hardware\+electrode\*.m');
fl(arrayfun(@(x) strcmp(x.name, 'probe.m'), fl) ) = [];
fl = arrayfun(@(x) x.name(1:end-2), fl, 'UniformOutput', false);
handles.popup_ChanMapLoad.Value  = 1;
handles.popup_ChanMapLoad.String = fl;

pr = handles.popup_ChanMapLoad.String{handles.popup_ChanMapLoad.Value};
handles.ops.probe = hardware.electrode.(pr);

% --- Default Headstage
fl = dir('+hardware\+headstage\*.m');
fl(arrayfun(@(x) strcmp(x.name, 'headstage.m'), fl) ) = [];
fl = arrayfun(@(x) x.name(1:end-2), fl, 'UniformOutput', false);

handles.popupHeadstageSelect.String = fl;
handles.popupHeadstageSelect.Value  = 1;

hs = handles.popupHeadstageSelect.String{handles.popupHeadstageSelect.Value};

handles.ops.headstage = hardware.headstage.(hs);

% --- Apply Channel map
if exist(handles.ops.chanMap, 'file')
    load(handles.ops.chanMap)
    handles.ops.chanMap = chanMap;
else
handles.ops.chanMap = handles.ops.headstage.mapChannels(handles.ops.probe);
end
% apply channel map
if numel(handles.ops.chanMap) < handles.ops.Nchan
    newChans = setdiff(1:handles.ops.Nchan, handles.ops.chanMap);
    handles.ops.chanMap = [handles.ops.chanMap(:); newChans(:)]';
end

% --- Default Referencing
handles.popup_ReferenceOptions.String = {'None', 'CommonAverage'};

handles.plot_buffer     = getBuffer(handles);
handles.plot_chanIdx    = 1:handles.ops.Nchan;

% --- Turn off panels until data is loaded
handles.Panel_ArtRem.Visible    = 'on';
handles.Panel_ChanMap.Visible   = 'on';
handles.Panel_Plot.Visible      = 'on';
handles.Panel_Reference.Visible = 'on';


for i = 1:handles.ops.Nchan
    handles.plot_handle(i) = plot(handles.plot_buffer(i,:), 'Color', repmat(.3, 1, 3)); hold on
end



% --- Scale bar and axes
handles.Ypos     = -100;
handles.axes1.YLim = [-100 handles.ops.Nchan*handles.plot_yoffset];
handles.scaleBarX = .9*diff(handles.axes1.XLim) + handles.axes1.XLim(1)*[1 1];
handles.scaleBarY = .05*diff(handles.axes1.YLim) + [0 200*handles.plot_yscale];
handles.scaleBarHandle = plot(handles.scaleBarX, handles.scaleBarY, 'k', 'Linewidth', 3);
handles.scaleBarText = text(handles.scaleBarX(1) + 50, mean(handles.scaleBarY), '200mV', 'Color', 'k');
grid(handles.axes1, 'on')
handles.axes1.XTick = '';
handles.axes1.YTick = '';
handles.axes1.Box = 'on';

% --- Slider steps
handles.Slider_Horizontal.Max = handles.fid_length;
handles.Slider_Horizontal.Value = handles.fid_position;
handles.Slider_Horizontal.SliderStep = [handles.plot_winWidth/handles.fid_length/100 handles.plot_winWidth/handles.fid_length/100];

guidata(hObject, handles);

handles = drawTraces(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in Button_xscaleminus.
function Button_xscaleminus_Callback(hObject, eventdata, handles)
% hObject    handle to Button_xscaleminus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = handles.plot_winWidth*1.1;
val = val - mod(val, 100);
handles.plot_winWidth   = floor(min(val, 2e5));
handles.plot_bufferSize = [handles.ops.Nchan handles.plot_winWidth];

fseek(handles.plot_fid, handles.fid_position, -1);

handles = drawTraces(handles);

handles.Slider_Horizontal.SliderStep = handles.ops.Nchan*[handles.plot_winWidth/handles.fid_length/100 handles.plot_winWidth/handles.fid_length/100];
% handles.Slider_Horizontal.SliderStep = [1 1]* (numel(handles.plot_buffer)/handles.fid_length); %/100); %handles.fid_length
guidata(hObject, handles)

% --- Executes on button press in Button_xscaleplus.
function Button_xscaleplus_Callback(hObject, ~, handles)
% hObject    handle to Button_xscaleplus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = handles.plot_winWidth*.9;
val = val - mod(val, 100);
handles.plot_winWidth   = floor(max(val, 100));
handles.plot_bufferSize = [handles.ops.Nchan handles.plot_winWidth];

fseek(handles.plot_fid, handles.fid_position, -1);

handles = drawTraces(handles);

handles.Slider_Horizontal.SliderStep = handles.ops.Nchan*[handles.plot_winWidth/handles.fid_length/100 handles.plot_winWidth/handles.fid_length/100];

guidata(hObject, handles)

% --- Executes on button press in Button_ElSepPlus.
function Button_ElSepPlus_Callback(hObject, eventdata, handles)
% hObject    handle to Button_ElSepPlus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.plot_yoffset = handles.plot_yoffset * 1.1;
fseek(handles.plot_fid, handles.fid_position, -1);
handles = drawTraces(handles);

guidata(hObject, handles);

% --- Executes on button press in Button_ElSepMinus.
function Button_ElSepMinus_Callback(hObject, eventdata, handles)
% hObject    handle to Button_ElSepMinus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.plot_yoffset = handles.plot_yoffset * .9;

fseek(handles.plot_fid, handles.fid_position, -1);
handles = drawTraces(handles);

guidata(hObject, handles);


% --- Executes on button press in Button_ZoomPlus.
function Button_ZoomPlus_Callback(hObject, eventdata, handles)
% hObject    handle to Button_ZoomPlus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.plot_yscale = handles.plot_yscale * 1.1;

fseek(handles.plot_fid, handles.fid_position, -1);

handles = drawTraces(handles);

guidata(hObject, handles);


% --- Executes on button press in Button_ZoomMinus.
function Button_ZoomMinus_Callback(hObject, eventdata, handles)
% hObject    handle to Button_ZoomMinus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.plot_yscale = handles.plot_yscale * .9;

fseek(handles.plot_fid, handles.fid_position, -1);

handles = drawTraces(handles);

guidata(hObject, handles);


% --- Executes on selection change in popup_ChanMapLoad.
function popup_ChanMapLoad_Callback(hObject, eventdata, handles)
% hObject    handle to popup_ChanMapLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_ChanMapLoad contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_ChanMapLoad
pr = handles.popup_ChanMapLoad.String{handles.popup_ChanMapLoad.Value};
handles.ops.probe   = hardware.electrode.(pr);
handles.ops.chanMap = handles.ops.headstage.mapChannels(handles.ops.probe);

fseek(handles.plot_fid, handles.fid_position, -1);
handles = drawTraces(handles);
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function popup_ChanMapLoad_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_ChanMapLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Input_ArtifactSize_Callback(hObject, eventdata, handles)
% hObject    handle to Input_ArtifactSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Input_ArtifactSize as text
%        str2double(get(hObject,'String')) returns contents of Input_ArtifactSize as a double
% val = str2double(get(hObject,'String'));
% handles.ops.artifactThresh = val/handles.ops.chHeaders(1).bitVolts;
% 
% % draw
% fseek(handles.plot_fid, handles.fid_position, -1);
% handles = drawTraces(handles);
% 
% % update GUI
% guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function Input_ArtifactSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Input_ArtifactSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Input_ArtifactChannelRange_Callback(hObject, eventdata, handles)
% hObject    handle to Input_ArtifactChannelRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Input_ArtifactChannelRange as text
%        str2double(get(hObject,'String')) returns contents of Input_ArtifactChannelRange as a double
% val = str2double(get(hObject,'String'));
% handles.ops.artifactNchans = val;
% 
% % draw
% fseek(handles.plot_fid, handles.fid_position, -1);
% handles = drawTraces(handles);
% 
% % update GUI
% guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function Input_ArtifactChannelRange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Input_ArtifactChannelRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_ReferenceOptions.
function popup_ReferenceOptions_Callback(hObject, eventdata, handles)
% hObject    handle to popup_ReferenceOptions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_ReferenceOptions contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_ReferenceOptions

% update Traces
fseek(handles.plot_fid, handles.fid_position, -1);
handles = drawTraces(handles);
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function popup_ReferenceOptions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_ReferenceOptions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Button_Close.
function Button_Close_Callback(hObject, eventdata, handles)
% hObject    handle to Button_Close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles, 'plot_fid')
    fclose(handles.plot_fid);
end

delete(hObject); % delete(hObject) closes the figure
close all force


% --- Executes on button press in checkbox_filter.
function checkbox_filter_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_filter
ops = handles.ops;

handles.filter.use =  get(hObject,'Value');
if handles.filter.use
    handles.ops.fshigh = str2double(handles.input_highcutoff.String);
    fslow = str2double(handles.input_lowpass.String);
    if fslow > 0 && fslow < ops.fs/2
        handles.ops.fslow = fslow;
    end
else
    handles.ops = rmfield(handles.ops, 'fshigh');
    if isfield(handles.ops, 'fslow')
        handles.ops = rmfield(handles.ops, 'fslow');
    end
end

handles = drawTraces(handles);
guidata(hObject, handles)


function input_lowpass_Callback(hObject, eventdata, handles)
% hObject    handle to input_lowpass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_lowpass as text
%        str2double(get(hObject,'String')) returns contents of input_lowpass as a double
% % val = str2double(get(hObject,'String'));
% % if val==0
% %     handles.ops = rmfield(handles.ops, 'fslow');
% % else
% %     handles.ops.fslow = val;
% % end
% % 
% % fseek(handles.plot_fid, handles.fid_position, -1);
% % handles = drawTraces(handles);
% % guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function input_lowpass_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_lowpass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function input_highcutoff_Callback(hObject, eventdata, handles)
% hObject    handle to input_highcutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_highcutoff as text
%        str2double(get(hObject,'String')) returns contents of input_highcutoff as a double
% val = str2double(get(hObject,'String'));
% 
% handles.ops.fshigh = val;
% fseek(handles.plot_fid, handles.fid_position, -1);
% handles = drawTraces(handles);
% guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function input_highcutoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_highcutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Button_Save.
function Button_Save_Callback(hObject, eventdata, handles)
% hObject    handle to Button_Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

chanMap   = handles.ops.chanMap;
n         = numel(handles.ops.chanMap);
connected = true(n, 1); % connected(1:2) = 0;
xcoords   = zeros(1,n);
ycoords   = 50*(1:n); 

% Often, multi-shank probes or tetrodes will be organized into groups of
% channels that cannot possibly share spikes with the rest of the probe. This helps
% the algorithm discard noisy templates shared across groups. In
% this case, we set kcoords to indicate which group the channel belongs to.
% In our case all channels are on the same shank in a single group so we
% assign them all to group 1. 

kcoords = ones(1,n);

% at this point in Kilosort we do data = data(connected, :), ycoords =
% ycoords(connected), xcoords = xcoords(connected) and kcoords =
% kcoords(connected) and no more channel map information is needed (in particular
% no "adjacency graphs" like in KlustaKwik). 
% Now we can save our channel map for the eMouse. 

% would be good to also save the sampling frequency here
fs = handles.ops.fs; 
ops = handles.ops;
ops.chanMap = fullfile(ops.root, 'chanMap.mat');
save(fullfile(handles.ops.root, 'ops.mat'), '-v7.3', '-struct', 'ops');
save(ops.chanMap, 'chanMap', 'connected', 'xcoords', 'ycoords', 'kcoords', 'fs')




% --- Executes on button press in Button_KiloSort.
function Button_KiloSort_Callback(hObject, eventdata, handles)
% hObject    handle to Button_KiloSort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fprintf('-----------------------------------------------------------------\n')
fprintf('-----------------------------------------------------------------\n')
fprintf('-----------------------------------------------------------------\n')
fprintf('Running Kilosort\n')

ops = handles.ops;
ops.parfor = true;

if ops.GPU
    fprintf('Using GPU for faster processing\n')
    fprintf('Opening gpuDevice. If this is the first time this was run during \nthis matlab session, this can be very slow\n')
    gpuDevice(1); % initialize GPU (will erase any existing GPU arrays)
end

[path_to, file, ext]=fileparts(ops.fbinary);
fOut=fullfile(path_to, [file '_hp' ext]);
% if ~exist(fOut, 'file')
    ops = removeArtifacts(ops);
% end
ops.fbinary = fOut;

% main spike-sorting routine
[rez, DATA, uproj] = preprocessData(ops); % preprocess data and extract spikes for initialization
rez                = fitTemplates(rez, DATA, uproj);  % fit templates iteratively
rez                = fullMPMU(rez, DATA);% extract final spike times (overlapping extraction)

% saving
fprintf('saving matlab results file\n')
save(fullfile(ops.root,  'rez.mat'), 'rez', 'ops', '-v7.3');

rez                = merge_posthoc2(rez);
fprintf('saving python files for Phy\n')

% save python results file for Phy
rezToPhy(rez, ops.root);

% remove temporary file
delete(ops.fproc);

% --- Executes on button press in checkbox_artifactremoval.
function checkbox_artifactremoval_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_artifactremoval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
useArtRemoval = get(hObject, 'Value');
if useArtRemoval
    val = str2double(handles.Input_ArtifactSize.String);
    handles.ops.artifactThresh = val/handles.ops.bitVolts;
    handles.ops.artifactNchans = str2double(handles.Input_ArtifactChannelRange.String);
else
    handles.ops.artifactThresh = inf;    
end
handles = drawTraces(handles);
guidata(hObject, handles)

% Hint: get(hObject,'Value') returns toggle state of checkbox_artifactremoval


% --- Executes on button press in Button_DetectWaveforms.
function Button_DetectWaveforms_Callback(hObject, eventdata, handles)
% hObject    handle to Button_DetectWaveforms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.figure1, 'renderer', 'painters')
saveas(handles.figure1, [datestr(now,'yyyymmdd_HHMM')], 'epsc')

% --- Executes on button press in Button_saveSpikeTimes.
function Button_saveSpikeTimes_Callback(hObject, eventdata, handles)
% hObject    handle to Button_saveSpikeTimes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% dat = zeros(handles.ops.Nchan,handles.fid_length/handles.ops.Nchan/2);
pos = 0;

fseek(handles.plot_fid, 0, 'bof');

sp=struct('st', [], 'ss', [], 'clu', []);

while ftell(handles.plot_fid) < handles.fid_length

    buffer = getBuffer(handles);
    s = {};
    id = {};
    for i = 1:handles.ops.Nchan
        buffer(i,:) = (buffer(i,:) - handles.plot_yoffset*i)/handles.plot_yscale;
        s{i} = detectSpikes(buffer(i,:), handles.ops.fs, -60);
        id{i} = ones(1,numel(s{i}))*i;
    end
    
    
    sp.ss  = [sp.ss cell2mat(s) + pos];
    sp.clu = [sp.clu cell2mat(id)];
    
    pos = pos + handles.plot_winWidth;
    
end

sp.st  = sp.ss/handles.ops.fs;

save(fullfile(handles.ops.root, 'sp.mat'), '-v7.3', '-struct', 'sp')





function [s, t] = detectSpikes(x, Fs, threshold)
% Detect spikes.
%   [s, t] = detectSpikes(x, Fs) detects spikes in x, where Fs the sampling
%   rate (in Hz). The outputs s and t are column vectors of spike times in
%   samples and ms, respectively. By convention the time of the zeroth
%   sample is 0 ms.

% detect local minima where at least one channel is above threshold
if nargin<3
    threshold = -100;
end
% noiseSD = median(abs(x)) / 0.6745;      % robust estimate of noise SD
% z = bsxfun(@rdivide, x, noiseSD);
% mz = min(z, [], 2);
% r = sqrt(sum(x .^ 2, 2));               % take norm for finding extrema
% dr = diff(r);
% s = find(mz(2 : end - 1) < threshold & dr(1 : end - 1) > 0 & dr(2 : end) < 0) + 1;
s = find(x(2 : end - 1) < threshold) + 1;
% s = s(s > 10 & s < size(x, 1) - 25);    % remove spikes close to boundaries

% % if multiple spikes occur within 1 ms we keep only the largest
% refractory = 1 / 1000 * Fs;
% N = numel(s);
% keep = true(N, 1);
% last = 1;
% for i = 2 : N
%     if s(i) - s(last) < refractory
%         if r(s(i)) > r(s(last))
%             keep(last) = false;
%             last = i;
%         else
%             keep(i) = false;
%         end
%     else
%         last = i;
%     end
% end
% s = s(keep);
t = s / Fs * 1000;                      % convert to real times in ms


% --- Executes on button press in checkbox_NotchFilter.
function checkbox_NotchFilter_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_NotchFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_NotchFilter
handles.ops.NotchFilter60 = get(hObject,'Value'); 
handles = drawTraces(handles);
guidata(hObject, handles)


% --- Executes on button press in button_MatlabSort.
function button_MatlabSort_Callback(hObject, eventdata, handles)
% hObject    handle to button_MatlabSort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupHeadstageSelect.
function popupHeadstageSelect_Callback(hObject, eventdata, handles)
% hObject    handle to popupHeadstageSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupHeadstageSelect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupHeadstageSelect


% --- Executes during object creation, after setting all properties.
function popupHeadstageSelect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupHeadstageSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

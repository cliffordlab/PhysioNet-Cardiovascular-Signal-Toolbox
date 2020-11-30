function varargout = ECGBeatFitter(varargin)
%
% ECGBeatFitter(ECG,Phase,ExpParamName),
% Graphical user interface for ECG approximation with Gaussian kernels.
%
% inputs:
% ECG: a single ECG waveform used for model training
% Phase: the phase corresponding to the ECG waveform
% ExpParamName (optional): The name of the vector containing the exported parameter (default:OptimumParams)
% ExpParamName (optional): The name of the vector containing the exported parameter (default:OptimumParams)
% Title (optional): The title of the plot
%
%
% Open Source ECG Toolbox, version 2.0, April 2008
% Released under the GNU General Public License
% Copyright (C) 2008  Reza Sameni
% Sharif University of Technology, Tehran, Iran -- GIPSA-LAB, INPG, Grenoble, France
% reza.sameni@gmail.com

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.

% ECGBEATFITTER M-file for ECGBeatFitter.fig
%      ECGBEATFITTER, by itself, creates a new ECGBEATFITTER or raises the existing
%      singleton*.
%
%      H = ECGBEATFITTER returns the handle to a new ECGBEATFITTER or the handle to
%      the existing singleton*.
%
%      ECGBEATFITTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ECGBEATFITTER.M with the given input arguments.
%
%      ECGBEATFITTER('Property','Value',...) creates a new ECGBEATFITTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ECGBeatFitter_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property
%      application
%      stop.  All inputs are passed to ECGBeatFitter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ECGBeatFitter

% Last Modified by GUIDE v2.5 08-Oct-2007 13:42:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ECGBeatFitter_OpeningFcn, ...
    'gui_OutputFcn',  @ECGBeatFitter_OutputFcn, ...
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


% --- Executes just before ECGBeatFitter is made visible.
function ECGBeatFitter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ECGBeatFitter (see VARARGIN)

%setappdata(hObject,'SelectedIndeces',0);
mn = varargin{1};
sd = varargin{2};
phase = varargin{3};
handles.mdlError = 'ModelErrorPercentage';
if (nargin>6),
    handles.ExpParamName = varargin{4};
else
    handles.ExpParamName = 'OptimumParams';
end

if (nargin>7),
    set(handles.text3,'String',varargin{5});
end

handles.Linehandles = [];
handles.Initmodhandle = [];
handles.Optmodhandle = [];

handles.SelectedIndeces = [];
handles.ECGmean = mn;
handles.ECGsd = sd;
handles.ECGphase = phase;
handles.OptimizedParameters = [];

axes(handles.axes1);

errorbar(mn,sd/2,'b');
hold on;
plot(mn,'r','linewidth',2);
grid on;
xlabel('Sample index');
ylabel('Amplitude (mV)');
set(handles.axes1,'ButtonDownFcn',@MyButtonDownFcn);
ECGcurve = get(handles.axes1,'Children');
set(ECGcurve,'ButtonDownFcn',@MyButtonDownFcn);
axis tight;
set(handles.figure1,'Pointer','crosshair');

% Choose default command line output for ECGBeatFitter
handles.output = hObject;
%handles.output = 0;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ECGBeatFitter wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = ECGBeatFitter_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

% varargout{1} = handles.output;
varargout{1} = handles.OptimizedParameters;

% --- Executes on button press in Reset.
function Reset_Callback(hObject, eventdata, handles)
% hObject    handle to Reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data = guidata(gcbo);

% delete(data.Initmodhandle);
delete(data.Optmodhandle);
% data.Initmodhandle = [];
data.Optmodhandle = [];

N = length(data.SelectedIndeces);
for i = 1:N
    delete(data.Linehandles(i));
end
data.SelectedIndeces = [];
data.Linehandles = [];

set(data.kernelnumber,'String','');
set(data.opterror,'String','');

guidata(gcbo,data);

% --- Executes on button press in Optimize.
function Optimize_Callback(hObject, eventdata, handles)
% hObject    handle to Optimize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data = guidata(gcbo);

% delete(data.Initmodhandle);
delete(data.Optmodhandle);

% localpeaks = zeros(size(data.ECGmean));
I = max(round(data.SelectedIndeces),1);
I = min(I,length(data.ECGmean));

P = length(I);
% localpeaks(I) = 1;

tetai = data.ECGphase(I(1:P));
alphai = 1.2*data.ECGmean(I(1:P));
bi = .04*ones(size(alphai));

options = optimset('TolX',1e-4,'TolFun',1e-4,'MaxIter',100);
InitParams = [alphai bi tetai];

OptParams = nlinfit(data.ECGphase,data.ECGmean,@ECGModel,InitParams,options);
% OptParams = lsqnonlin(@(InitParams) ECGModelError(InitParams,ECGmn,Phasemn,0),InitParams,InitParams-2,InitParams+2,options);
% Model0 = ECGModelError(InitParams,data.ECGmean,data.ECGphase,1);
Model = ECGModelError(OptParams,data.ECGmean,data.ECGphase,1);

ax = axis;
% data.Initmodhandle = plot(Model0,'g');
data.Optmodhandle = plot(Model,'k','LineWidth',2);

L = length(OptParams);

ph = OptParams(2*L/3+1:L);
for i = 1:L/3,
    if(ph(i)<min(data.ECGphase) || ph(i)>max(data.ECGphase))
        ph(i) = mean(data.ECGphase); % place the possible erroneous estimate at the middle of the ECG beat (close to the R-peak)
    end

    pos = find(data.ECGphase>=ph(i),1);

    I = max(pos,1);
    I = min(I,length(data.ECGphase));
    data.SelectedIndeces(i) = I;
    set(data.Linehandles(i),'XData',[I I]);
end
% % %
% % % P = length(I);
% % % % localpeaks(I) = 1;
% % %
% % % tetai = data.ECGphase(I(1:P));



axis(ax);

data.OptimizedParameters = OptParams;

er = 100*mean((Model-data.ECGmean).^2)/mean(data.ECGmean.^2);

set(data.kernelnumber,'String',length(data.SelectedIndeces));
set(data.opterror,'String',round(100*er)/100);

% Modified 25/1/2019 by Reza Sameni:
% assignin('base',handles.mdlError,er);
assignin('caller',handles.mdlError,er);


guidata(gcbo,data);

% --- Executes on button press in ExportData.
function ExportData_Callback(hObject, eventdata, handles)
% hObject    handle to ExportData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Modified 25/1/2019 by Reza Sameni:
% assignin('base',handles.ExpParamName,handles.OptimizedParameters);
assignin('caller',handles.ExpParamName,handles.OptimizedParameters);

% handles.output = 1; % Function has finished execution
% guidata(hObject, handles);
uiresume(handles.figure1);


function MyButtonDownFcn(hObject, eventdata, handles)
pn = get(gca,'Currentpoint');
% hndl = plot(pn(1,1),pn(1,2),'ro');
ax = axis;
ln = line([pn(1,1),pn(1,1)],[ax(3),ax(4)]);
set(ln,'Color',[0 1 .5]);

data = guidata(gcbo);
data.SelectedIndeces = [data.SelectedIndeces pn(1,1)];
data.Linehandles = [data.Linehandles ln];

guidata(gcbo,data);


% --------------------------------------------------------------------
function Helpmenu_Callback(hObject, eventdata, handles)
% hObject    handle to Helpmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function About_1_Callback(hObject, eventdata, handles)
% hObject    handle to About_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
AboutPath = which('About.htm');
web(AboutPath);

% --------------------------------------------------------------------
function Help_1_Callback(hObject, eventdata, handles)
% hObject    handle to Help_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

HelpPath = which('Help.htm');
web(HelpPath);


% --- Executes on button press in Exit.
function Exit_Callback(hObject, eventdata, handles)
% hObject    handle to Exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% assignin('base',handles.ExpParamName,[]);
uiresume(handles.figure1);

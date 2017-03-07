function varargout = blood_flow(varargin)
% BLOOD_FLOW M-file for blood_flow.fig
%      BLOOD_FLOW, by itself, creates a new BLOOD_FLOW or raises the existing
%      singleton*.
%
%      H = BLOOD_FLOW returns the handle to a new BLOOD_FLOW or the handle to
%      the existing singleton*.
%
%      BLOOD_FLOW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BLOOD_FLOW.M with the given input arguments.
%
%      BLOOD_FLOW('Property','Value',...) creates a new BLOOD_FLOW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before blood_flow_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to blood_flow_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help blood_flow

% Last Modified by GUIDE v2.5 09-Nov-2015 15:33:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @blood_flow_OpeningFcn, ...
                   'gui_OutputFcn',  @blood_flow_OutputFcn, ...
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


% --- Executes just before blood_flow is made visible.
function blood_flow_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to blood_flow (see VARARGIN)

% Choose default command line output for blood_flow
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes blood_flow wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = blood_flow_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
aL=str2double(get(handles.edit3,'String'));
aR=str2double(get(handles.edit4,'String'));
uL=str2double(get(handles.edit1,'String'));
uR=str2double(get(handles.edit2,'String'));
NCELLS=str2double(get(handles.edit5,'String'));
CFL=str2double(get(handles.edit6,'String'));
blood_flow_compute(aL,aR,uL,uR,NCELLS,CFL,handles)

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.figure1)


function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
aL=str2double(get(handles.edit3,'String'));
aR=str2double(get(handles.edit4,'String'));
uL=str2double(get(handles.edit1,'String'));
uR=str2double(get(handles.edit2,'String'));
NCELLS=str2double(get(handles.edit5,'String'));
CFL=str2double(get(handles.edit6,'String'));

if(aL>0 && aR>0 && ~isnan(uL) && ~isnan(uR) && NCELLS>0 && CFL>0)
    set(handles.pushbutton1,'Enable','on')
else
    set(handles.pushbutton1,'Enable','off')
end

% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double
aL=str2double(get(handles.edit3,'String'));
aR=str2double(get(handles.edit4,'String'));
uL=str2double(get(handles.edit1,'String'));
uR=str2double(get(handles.edit2,'String'));
NCELLS=str2double(get(handles.edit5,'String'));
CFL=str2double(get(handles.edit6,'String'));

if(aL>0 && aR>0 && ~isnan(uL) && ~isnan(uR) && NCELLS>0 && CFL>0)
    set(handles.pushbutton1,'Enable','on')
else
    set(handles.pushbutton1,'Enable','off')
end

% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
aL=str2double(get(handles.edit3,'String'));
aR=str2double(get(handles.edit4,'String'));
uL=str2double(get(handles.edit1,'String'));
uR=str2double(get(handles.edit2,'String'));
NCELLS=str2double(get(handles.edit5,'String'));
CFL=str2double(get(handles.edit6,'String'));

if(aL>0 && aR>0 && ~isnan(uL) && ~isnan(uR) && NCELLS>0 && CFL>0)
    set(handles.pushbutton1,'Enable','on')
else
    set(handles.pushbutton1,'Enable','off')
end

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
aL=str2double(get(handles.edit3,'String'));
aR=str2double(get(handles.edit4,'String'));
uL=str2double(get(handles.edit1,'String'));
uR=str2double(get(handles.edit2,'String'));
NCELLS=str2double(get(handles.edit5,'String'));
CFL=str2double(get(handles.edit6,'String'));

if(aL>0 && aR>0 && ~isnan(uL) && ~isnan(uR) && NCELLS>0 && CFL>0)
    set(handles.pushbutton1,'Enable','on')
else
    set(handles.pushbutton1,'Enable','off')
end

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
aL=str2double(get(handles.edit3,'String'));
aR=str2double(get(handles.edit4,'String'));
uL=str2double(get(handles.edit1,'String'));
uR=str2double(get(handles.edit2,'String'));
NCELLS=str2double(get(handles.edit5,'String'));
CFL=str2double(get(handles.edit6,'String'));

if(aL>0 && aR>0 && ~isnan(uL) && ~isnan(uR) && NCELLS>0 && CFL>0)
    set(handles.pushbutton1,'Enable','on')
else
    set(handles.pushbutton1,'Enable','off')
end

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double
aL=str2double(get(handles.edit3,'String'));
aR=str2double(get(handles.edit4,'String'));
uL=str2double(get(handles.edit1,'String'));
uR=str2double(get(handles.edit2,'String'));
NCELLS=str2double(get(handles.edit5,'String'));
CFL=str2double(get(handles.edit6,'String'));

if(aL>0 && aR>0 && ~isnan(uL) && ~isnan(uR) && NCELLS>0 && CFL>0)
    set(handles.pushbutton1,'Enable','on')
else
    set(handles.pushbutton1,'Enable','off')
end

% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes during object creation, after setting all properties.
function axes9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes9

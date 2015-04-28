function varargout = human_recognition_test(varargin)

% HUMAN_RECOGNITION_TEST MATLAB code for human_recognition_test.fig
%      HUMAN_RECOGNITION_TEST, by itself, creates a new HUMAN_RECOGNITION_TEST or raises the existing
%      singleton*.
%
%      H = HUMAN_RECOGNITION_TEST returns the handle to a new HUMAN_RECOGNITION_TEST or the handle to
%      the existing singleton*.
%
%      HUMAN_RECOGNITION_TEST('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HUMAN_RECOGNITION_TEST.M with the given input arguments.
%
%      HUMAN_RECOGNITION_TEST('Property','Value',...) creates a new HUMAN_RECOGNITION_TEST or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before human_recognition_test_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to human_recognition_test_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help human_recognition_test

% Last Modified by GUIDE v2.5 27-Apr-2015 23:34:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @human_recognition_test_OpeningFcn, ...
                   'gui_OutputFcn',  @human_recognition_test_OutputFcn, ...
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


% --- Executes just before human_recognition_test is made visible.
function human_recognition_test_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to human_recognition_test (see VARARGIN)

% Choose default command line output for human_recognition_test
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


image_path = 'training\';
image_list =  dir([image_path '*jpg']);
n_images = size(image_list, 1);

warped_image_path = 'warped_images\';
warp_list = dir([warped_image_path '*jpg']);


global image_data;
image_data = cell(n_images);
global warped_image_data;
warped_image_data = cell(n_images);

global counter;
counter = 1;

for j=1:n_images
       j
       filename = strcat(image_path, image_list(j).name);
       X = double(imread(filename))./255;
       image_data{j} = X;
end

for j=1:n_images
       j
       filename = strcat(warped_image_path, warp_list(j).name);
       X = double(imread(filename))./255;
       warped_image_data{j} = X;
end

 axes(handles.axes1);
 imshow(image_data{counter});
 
 axes(handles.axes8);
 imshow(warped_image_data{counter});
 
 set([handles.axes2, handles.axes3,handles.axes4, handles.axes5,handles.axes6, handles.axes7], 'buttondownfcn', @choose_image);
 
 
 
 function choose_image(hObject, eventdata)
  global warped_image_data;
  global image_data;
  global counter;
  counter = counter +1;
  handles = guidata(hObject);
  theaxes = ancestor(hObject,'axes');
  if theaxes == handles.axes2
    %axes2 stuff
     
  elseif theaxes == handles.axes3
    %axes3 stuff
  elseif theaxes == handles.axes4
    %axes4 stuff
  elseif theaxes == handles.axes5
    %axes5 stuff
  elseif theaxes == handles.axes6
    %axes6 stuff
  elseif theaxes == handles.axes7
    %axes7 stuff
  end
  axes(handles.axes1);
  imshow(image_data{counter});

  axes(handles.axes8);
  imshow(warped_image_data{counter});
  set([handles.axes2, handles.axes3,handles.axes4, handles.axes5,handles.axes6, handles.axes7], 'buttondownfcn', @choose_image);
     
     
        
 
 



% UIWAIT makes human_recognition_test wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = human_recognition_test_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2


% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4


% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton5


% --- Executes on button press in radiobutton6.
function radiobutton6_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton6


% --- Executes on button press in radiobutton7.
function radiobutton7_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton7


% --- Executes on button press in radiobutton8.
function radiobutton8_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton8


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



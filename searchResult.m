function varargout = searchResult(varargin)
% SEARCHRESULT MATLAB code for searchResult.fig
%      SEARCHRESULT, by itself, creates a new SEARCHRESULT or raises the existing
%      singleton*.
%
%      H = SEARCHRESULT returns the handle to a new SEARCHRESULT or the handle to
%      the existing singleton*.
%
%      SEARCHRESULT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEARCHRESULT.M with the given input arguments.
%
%      SEARCHRESULT('Property','Value',...) creates a new SEARCHRESULT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before searchResult_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to searchResult_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help searchResult

% Last Modified by GUIDE v2.5 11-Dec-2017 17:17:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @searchResult_OpeningFcn, ...
                   'gui_OutputFcn',  @searchResult_OutputFcn, ...
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


% --- Executes just before searchResult is made visible.
function searchResult_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to searchResult (see VARARGIN)

% Choose default command line output for searchResult
handles.output = hObject;

h = findobj('Tag','Gui1');

if ~isempty(h)
    
    files = getappdata(h,'files');
    ids = getappdata(h,'ids');
    I = getappdata(h,'qim');
    imtitle = getappdata(h,'imtitle');
    resultACC = getappdata(h,'resultACC');
       
    axes(handles.axes2);
    imshow(I);
    title(imtitle,'Interpreter','none');
    
    set( handles.text4,'String',sprintf('Accuracy: %f', resultACC) );
    drawnow;
    
    axes(handles.axes1);
    fid = fopen('oxford\groundtruth\rank_list.txt', 'w');    
    for i=1:size(ids{1},2)
        % Show only 10 highest score images
        if i<=10            
            subplot(3, 5, 5+i); 
            imshow(imread(fullfile('oxford\images\', files(ids{1}(i)).name)));
            title(files(ids{1}(i)).name,'Interpreter','none');
        end
        fprintf(fid, '%s\n', files(ids{1}(i)).name(1:end-4));
    end
    fclose(fid);
    
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes searchResult wait for user response (see UIRESUME)
% uiwait(handles.Gui2);


% --- Outputs from this function are returned to the command line.
function varargout = searchResult_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout{1} = handles.output;

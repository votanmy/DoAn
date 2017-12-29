function varargout = main(varargin)
% MAIN MATLAB code for main.fig
%      MAIN, by itself, creates a new MAIN or raises the existing
%      singleton*.
%
%      H = MAIN returns the handle to a new MAIN or the handle to
%      the existing singleton*.
%
%      MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN.M with the given input arguments.
%
%      MAIN('Property','Value',...) creates a new MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before main_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help main

% Last Modified by GUIDE v2.5 28-Dec-2017 19:39:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @main_OpeningFcn, ...
                   'gui_OutputFcn',  @main_OutputFcn, ...
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


% --- Executes just before main is made visible.
function main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to main (see VARARGIN)

% Choose default command line output for main
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes main wait for user response (see UIRESUME)
% uiwait(handles.Gui0);


% --- Outputs from this function are returned to the command line.
function varargout = main_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

pause(0.001);

%% init parameter
global if_weight if_norm if_dist files
addpath('AKM');
run('vlfeat\toolbox\vl_setup.m');
datasetDir = 'oxford\images\';
num_words = 1000;
num_iterations = 5;
num_trees = 8;
%dim = 128;
if_weight = 'tfidf';
if_norm = 'l1';
if_dist = 'l1';
%verbose=1;
%% Compute SIFT features
if ~exist('oxford\feat\feature.bin', 'file')
    %fprintf('Computing SIFT features:\n');
    set(handles.text4,'String','Computing SIFT features......');    
    drawnow;
    features = zeros(128, 2000000);
    nfeat = 0;
    files = dir(fullfile(datasetDir, '*.jpg'));
    nfiles = length(files);
    features_per_image = zeros(1,nfiles);
    for i=1:nfiles
        %fprintf('Extracting features %d/%d\n', i, nfiles);
        set(handles.text4,'String','Extracting features......');
        drawnow;
        imgPath = strcat(datasetDir, files(i).name);
        I = im2single(rgb2gray(imread(imgPath)));
        I = imresize(I, 0.6);
        [frame, sift] = vl_covdet(I, 'method', 'Hessian', 'estimateAffineShape', true);
        
        if nfeat+size(sift,2) > size(features,2)
            features = [features zeros(128,1000000)];
        end
        features(:,nfeat+1:nfeat+size(sift,2)) = sift;
        nfeat = nfeat+size(sift,2);
        features_per_image(i) = size(sift, 2);
    end
    features = features(:,1:nfeat);
    fid = fopen('oxford\feat\feature.bin', 'w');
    fwrite(fid, features, 'float64');
    fclose(fid);
    
    save('oxford\feat\feat_info.mat', 'features_per_image', 'files');
else
    %fprintf('Loading SIFT features:\n');    
    set(handles.text4,'String','Loading SIFT features......');
    drawnow;
    file = dir('oxford\feat\feature.bin');
    %features = zeros(128, file.bytes/(4*128), 'single');

    fid = fopen('oxford\feat\feature.bin', 'r');
    features = fread(fid, [128, file.bytes/(4*128)], 'float64');
    fclose(fid);
    
    load('oxford\feat\feat_info.mat');
end

%% compute rootSIFT
% %fprintf('Computing rootSIFT features:\n');
% set(handles.text4,'String','Computing rootSIFT features......');
% drawnow;
% num_features = size(features, 2);
% %rootSIFT = zeros(dim, num_features);
% 
% poolObj = parpool('local', 4);
% for k = 1:5000000:num_features
%     eIdx = k+5000000-1;
%     if eIdx > num_features
%         eIdx = num_features;
%     end
%     parfor i=k:eIdx
%         features(:, i) = sqrt(features(:, i) / sum(features(:,i)));
%     end
% end
% delete(poolObj);

%% Run AKM to build dictionary
global dict_words inv_file dict
%fprintf('Building the dictionary:\n');
set(handles.text5,'String','Building the dictionary......'); 
drawnow;
num_images = length(files);
dict_params =  {num_iterations, 'kdt', num_trees};

% build the dictionary
if exist('oxford\feat\dict.mat', 'file')
    load('oxford\feat\dict.mat');
else
    randIndex = randperm(size(features,2));
    dict_words = ccvBowGetDict(features(:,randIndex(1:100000)), [], [], num_words, 'flat', 'akmeans',...
        [], dict_params);    % ccvRanSeek.m in AKM folder has been changed on line 26
    save('oxford\feat\dict.mat', 'dict_words');
end

% compute sparse frequency vector
%fprintf('Computing the words\n');
set(handles.text6,'String','Computing the words......');  
drawnow;
dict = ccvBowGetWordsInit(dict_words, 'flat', 'akmeans', [], dict_params);

if exist('oxford\feat\words.mat', 'file')
    load('oxford\feat\words.mat');
else
    words = cell(1, num_images);
    for i=1:num_images
        %fprintf('Quantizing %d/%d images\n', i, num_images);
        set(handles.text6,'String','Quantizing......'); 
        drawnow;
        if i==1
            bIndex = 1;
        else
            bIndex = sum(features_per_image(1:i-1))+1;
        end
        eIndex = bIndex + features_per_image(i)-1;
        words{i} = ccvBowGetWords(dict_words, features(:, bIndex:eIndex), [], dict);
    end;
    save('oxford\feat\words.mat', 'words');
end
% fprintf('Computing sparse frequency vector\n');
% dict = ccvBowGetWordsInit(dict_words, 'flat', 'akmeans', [], dict_params);
% words = ccvBowGetWords(dict_words, root_sift, [], dict);
% ccvBowGetWordsClean(dict);

% create an inverted file for the images
%fprintf('Creating and searching an inverted file\n');
set(handles.text7,'String','Creating and searching an inverted file......');    
drawnow;
inv_file = ccvInvFileInsert([], words, num_words);
ccvInvFileCompStats(inv_file, if_weight, if_norm);
%save('inverted_file.mat', 'if_weight', 'if_norm', 'if_dist', 'inv_file');

set(handles.text8,'String','Stop loading. Start searching....');    
drawnow;
guidata(hObject,handles);
setappdata(handles.Gui0,'files',files);
setappdata(handles.Gui0,'dict_words',dict_words);
setappdata(handles.Gui0,'dict',dict); 
setappdata(handles.Gui0,'inv_file',inv_file);    
setappdata(handles.Gui0,'if_weight',if_weight);
setappdata(handles.Gui0,'if_norm',if_norm);
setappdata(handles.Gui0,'if_dist',if_dist);
run('imsearch.m');

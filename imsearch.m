function varargout = imsearch(varargin)
% IMSEARCH MATLAB code for imsearch.fig
%      IMSEARCH, by itself, creates a new IMSEARCH or raises the existing
%      singleton*.
%
%      H = IMSEARCH returns the handle to a new IMSEARCH or the handle to
%      the existing singleton*.
%
%      IMSEARCH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMSEARCH.M with the given input arguments.
%
%      IMSEARCH('Property','Value',...) creates a new IMSEARCH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before imsearch_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to imsearch_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help imsearch

% Last Modified by GUIDE v2.5 29-Dec-2017 06:12:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @imsearch_OpeningFcn, ...
                   'gui_OutputFcn',  @imsearch_OutputFcn, ...
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


% --- Executes just before imsearch is made visible.
function imsearch_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to imsearch (see VARARGIN)

% Choose default command line output for imsearch
handles.output = hObject;



% %% init parameter
% global if_weight if_norm if_dist files
% addpath('AKM');
% run('vlfeat\toolbox\vl_setup.m');
% datasetDir = 'oxford\images\';
% num_words = 1000;
% num_iterations = 5;
% num_trees = 8;
% dim = 128;
% if_weight = 'tfidf';
% if_norm = 'l1';
% if_dist = 'l1';
% verbose=1;
% %% Compute SIFT features
% if ~exist('oxford\feat\feature.bin', 'file')
%     fprintf('Computing SIFT features:\n');
%     
%     features = zeros(128, 2000000);
%     nfeat = 0;
%     files = dir(fullfile(datasetDir, '*.jpg'));
%     nfiles = length(files);
%     features_per_image = zeros(1,nfiles);
%     for i=1:nfiles
%         fprintf('Extracting features %d/%d\n', i, nfiles);
%         imgPath = strcat(datasetDir, files(i).name);
%         I = im2single(rgb2gray(imread(imgPath)));
%         I = imresize(I, 0.6);
%         [frame, sift] = vl_covdet(I, 'method', 'Hessian', 'estimateAffineShape', true);
%         
%         if nfeat+size(sift,2) > size(features,2)
%             features = [features zeros(128,1000000)];
%         end
%         features(:,nfeat+1:nfeat+size(sift,2)) = sift;
%         nfeat = nfeat+size(sift,2);
%         features_per_image(i) = size(sift, 2);
%     end
%     features = features(:,1:nfeat);
%     fid = fopen('oxford\feat\feature.bin', 'w');
%     fwrite(fid, features, 'float64');
%     fclose(fid);
%     
%     save('oxford\feat\feat_info.mat', 'features_per_image', 'files');
% else
%     fprintf('Loading SIFT features:\n');
%     file = dir('oxford\feat\feature.bin');
%     %features = zeros(128, file.bytes/(4*128), 'single');
% 
%     fid = fopen('oxford\feat\feature.bin', 'r');
%     features = fread(fid, [128, file.bytes/(4*128)], 'float64');
%     fclose(fid);
%     
%     load('oxford\feat\feat_info.mat');
% end
% 
% %% compute rootSIFT
% fprintf('Computing rootSIFT features:\n');
% num_features = size(features, 2);
% %rootSIFT = zeros(dim, num_features);
% 
% % poolObj = parpool('local', 4);
% % for k = 1:5000000:num_features
% %     eIdx = k+5000000-1;
% %     if eIdx > num_features
% %         eIdx = num_features;
% %     end
% %     parfor i=k:eIdx
% %         features(:, i) = sqrt(features(:, i) / sum(features(:,i)));
% %     end
% % end
% % delete(poolObj);
% 
% %% Run AKM to build dictionary
% global dict_words inv_file dict
% fprintf('Building the dictionary:\n');
% num_images = length(files);
% dict_params =  {num_iterations, 'kdt', num_trees};
% 
% % build the dictionary
% if exist('oxford\feat\dict.mat', 'file')
%     load('oxford\feat\dict.mat');
% else
%     randIndex = randperm(size(features,2));
%     dict_words = ccvBowGetDict(features(:,randIndex(1:100000)), [], [], num_words, 'flat', 'akmeans',...
%         [], dict_params);    % ccvRanSeek.m in AKM folder has been changed on line 26
%     save('oxford\feat\dict.mat', 'dict_words');
% end
% 
% % compute sparse frequency vector
% fprintf('Computing the words\n');
% dict = ccvBowGetWordsInit(dict_words, 'flat', 'akmeans', [], dict_params);
% 
% if exist('oxford\feat\words.mat', 'file')
%     load('oxford\feat\words.mat');
% else
%     words = cell(1, num_images);
%     for i=1:num_images
%         fprintf('Quantizing %d/%d images\n', i, num_images);
%         if i==1
%             bIndex = 1;
%         else
%             bIndex = sum(features_per_image(1:i-1))+1;
%         end
%         eIndex = bIndex + features_per_image(i)-1;
%         words{i} = ccvBowGetWords(dict_words, features(:, bIndex:eIndex), [], dict);
%     end;
%     save('oxford\feat\words.mat', 'words');
% end
% %fprintf('Computing sparse frequency vector\n');
% %dict = ccvBowGetWordsInit(dict_words, 'flat', 'akmeans', [], dict_params);
% %words = ccvBowGetWords(dict_words, root_sift, [], dict);
% %ccvBowGetWordsClean(dict);
% 
% % create an inverted file for the images
% fprintf('Creating and searching an inverted file\n');
% inv_file = ccvInvFileInsert([], words, num_words);
% ccvInvFileCompStats(inv_file, if_weight, if_norm);
% %save('inverted_file.mat', 'if_weight', 'if_norm', 'if_dist', 'inv_file');




% Update handles structure
guidata(hObject, handles);

% UIWAIT makes imsearch wait for user response (see UIRESUME)
% uiwait(handles.Gui1);


% --- Outputs from this function are returned to the command line.
function varargout = imsearch_OutputFcn(hObject, eventdata, handles) 
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

global slectedImg Icr

if exist('cropim.jpg', 'file')
    delete 'cropim.jpg';
end
[imgname,~] = imgetfile('InitialPath','oxford/query_images/');
[~,slectedImg,~] = fileparts(imgname);
im = imread(imgname);
im = im2double(im); %converts to double

axes(handles.axes1);
imshow(im);
set(handles.text4,'String',slectedImg);

Icr = imcrop(im);
pause(2);
if ~isempty(Icr)
    imshow(Icr, []);
    set(handles.text4,'String',strcat(slectedImg, ' Cropped'));
    imwrite(Icr,'cropim.jpg');
end

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h = findobj('Tag','Gui0');
files = getappdata(h,'files');
dict_words = getappdata(h,'dict_words');
dict = getappdata(h,'dict');
inv_file = getappdata(h,'inv_file');
if_weight = getappdata(h,'if_weight');
if_norm = getappdata(h,'if_norm');
if_dist = getappdata(h,'if_dist');

%% Query images
global slectedImg
q_files = dir(fullfile('oxford\groundtruth', '*query.txt'));
%oxc1_all_souls_000013 136.5 34.1 648.5 955.7
ntop = 0;

% load query image
for k=1:length(q_files)

    fid = fopen(strcat('oxford\groundtruth\', q_files(k).name), 'r');
    str = fgetl(fid);
    [image_name, remain] = strtok(str, ' ');
    fclose(fid);
    numbers = str2num(remain);
    if strcmp( slectedImg, image_name(6:end) )
        
        if exist('cropim.jpg', 'file')
            Icr = imread('cropim.jpg');
            if size(Icr,3)==3
                I = im2single(rgb2gray(Icr));
            else
                I = im2single(Icr);
            end
            
            %imgIF = imfinfo(Icr);
            [imgW,imgH] = size(I);
            x1 = 0;
            y1 = 0;
            x2 = imgW;
            y2 = imgH;
            
            imtitle = strcat(slectedImg, ' Cropped');
            % compute rootSIFT features
            [frame, sift] = vl_covdet(I, 'method', 'Hessian', 'estimateAffineShape', true); 
            sift = sift(:,(frame(1,:)<=x2) &  (frame(1,:) >= x1) & (frame(2,:) <= y2) & (frame(2,:) >= y1));
            
            % Test on an image
            %global files dict_words dict inv_file if_weight if_norm if_dist
            q_words = cell(1,1);
            q_words{1} = ccvBowGetWords(dict_words, double(sift), [], dict);
            [ids dists] = ccvInvFileSearch(inv_file, q_words(1), if_weight, if_norm, if_dist, ntop);
            
            fid = fopen('oxford\groundtruth\rank_list_cropped.txt', 'w');
            for i=1:size(ids{1},2)
                % Build the rank list for the cropped image
                fprintf(fid, '%s\n', files(ids{1}(i)).name(1:end-4));
            end
            fclose(fid);

            script = ['oxford\groundtruth\Test.exe oxford\groundtruth\', ...
                q_files(k).name(1:end-10), ...
                ' oxford\groundtruth\rank_list_cropped.txt',...
                ' >oxford\cropped_result\', image_name(6:end), '_crop_result.txt']; %q_files(k).name(1:end-10)
            system(script);
            
            result_file = fopen(strcat('oxford\cropped_result\', image_name(6:end), '_crop_result.txt'),'r');
        
        else 
            x1 = numbers(1);
            y1 = numbers(2);
            x2 = numbers(3);
            y2 = numbers(4);
            file = strcat('oxford\images\', image_name(6:end), '.jpg');  
            
            I = im2single(rgb2gray(imread(file)));
            imtitle = image_name(6:end);
            % compute rootSIFT features
            %[~, sift] = vl_covdet(I, 'method', 'Hessian', 'estimateAffineShape', true);
            [frame, sift] = vl_covdet(I, 'method', 'Hessian', 'estimateAffineShape', true);
            sift = sift(:,(frame(1,:)<=x2) &  (frame(1,:) >= x1) & (frame(2,:) <= y2) & (frame(2,:) >= y1));
            
            % Test on an image
            %global files dict_words dict inv_file if_weight if_norm if_dist
            q_words = cell(1,1);
            q_words{1} = ccvBowGetWords(dict_words, double(sift), [], dict);
            [ids dists] = ccvInvFileSearch(inv_file, q_words(1), if_weight, if_norm, if_dist, ntop);
            
            fid = fopen('oxford\groundtruth\rank_list.txt', 'w');
            for i=1:size(ids{1},2)
                % Build the rank list slected image
                fprintf(fid, '%s\n', files(ids{1}(i)).name(1:end-4));
            end
            fclose(fid);            

            script = ['oxford\groundtruth\Test.exe oxford\groundtruth\', ...
                q_files(k).name(1:end-10), ...
                ' oxford\groundtruth\rank_list.txt',...
                ' >oxford\result\', image_name(6:end), '_result.txt']; %q_files(k).name(1:end-10)
            system(script);
            
            result_file = fopen(strcat('oxford\result\', image_name(6:end), '_result.txt'),'r');
            
        end
        
        resultACC = fscanf(result_file, '%f');

        %imshow(I); hold on;
        %plot([x1 x2], [y1 y1], 'g');
        %plot([x2 x2], [y1 y2], 'g');
        %plot([x2 x1], [y2 y2], 'g');
        %plot([x1 x1], [y2 y1], 'g');
        %hold off;
        % compute rootSIFT features
        %[frame, sift] = vl_covdet(I, 'method', 'Hessian', 'estimateAffineShape', true);
        %sift = sift(:,(frame(1,:)<=x2) &  (frame(1,:) >= x1) & (frame(2,:) <= y2) & (frame(2,:) >= y1)); 

        % Test on an image
%         global files dict_words dict inv_file if_weight if_norm if_dist
%         q_words = cell(1,1);
%         q_words{1} = ccvBowGetWords(dict_words, double(sift), [], dict);
%         [ids dists] = ccvInvFileSearch(inv_file, q_words(1), if_weight, if_norm, if_dist, ntop);
% 
%         script = ['oxford\groundtruth\Test.exe oxford\groundtruth\', ...
%             q_files(k).name(1:end-10), ...
%             ' oxford\groundtruth\rank_list.txt',...
%             ' >oxford\result\', image_name(6:end), '_result.txt']; %q_files(k).name(1:end-10)
%         system(script);         
        
    else
        continue;
    end
end
   
    
guidata(hObject,handles);
setappdata(handles.Gui1,'files',files);
setappdata(handles.Gui1,'ids',ids);
setappdata(handles.Gui1,'imtitle',imtitle); 
setappdata(handles.Gui1,'qim',I); 
setappdata(handles.Gui1,'resultACC',resultACC); 
  
run('searchResult.m');

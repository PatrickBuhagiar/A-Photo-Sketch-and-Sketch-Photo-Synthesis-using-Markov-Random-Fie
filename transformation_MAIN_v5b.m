%% Clear variables

clearvars;


%% Initialise

CLK = clock;
fprintf('\nFunction started at %02d/%02d/%d, %d:%02.0f:%02.0f \n', CLK(3),CLK(2),CLK(1),CLK(4),CLK(5),CLK(6) );
tStart = tic;





%% 0.1 Add the required paths
% addpath('F:\Google Drive\MATLAB\eig_trans__Reuben\MATLAB\');
% addpath('F:\Google Drive\MATLAB\eig_trans__Reuben\MATLAB\dataset_parsing');
% addpath('F:\Google Drive\MATLAB\eig_trans__Reuben\MATLAB\general');
% addpath('F:\Google Drive\MATLAB\eig_trans__Reuben\MATLAB\registration');
%
%addpath('F:\Google Drive\MATLAB\eig_trans__Reuben\MATLAB\eigentransformaiton');

%base_path = 'F:\Google Drive\MATLAB\PhD\My Work\Transformations and Face Matchers\';

%addFiles

%% 0.2 Get information from files or load MAT file containing data: normalised photo and sketch images etc

% Database dependent;

% This .mat file must have 'data_photos_sketches', which must contain 'photos_norm',
% 'sketches_norm' and corresponding information (subject ID, eye
% coordinates etc)

param.DB_name = 'ColorFERET';

if strcmpi(param.DB_name, 'ColorFERET')
    %mat_filename = 'F:\Google Drive\MATLAB\Eigentransformation\ColorFERET_norm_ALL_and_sketches__renamed2_and_using_CUFSF_file_names2_newNamingConvention2.mat';
    mp = mfilename('fullpath');
    [pathstr, name, ext] = fileparts(mp);
    %mat_filename = 'F:\Google Drive\MATLAB\PhD\My Work\Transformations and Face Matchers\Eigentransformation and Eigenfaces\ColorFERET_norm_ALL_and_sketches__renamed2_and_using_CUFSF_file_names2_newNamingConvention2.mat'; %([pathstr '\..\' 'Data\ColorFERET_norm_ALL_and_sketches__renamed2_and_using_CUFSF_file_names2_newNamingConvention2.mat']);%'F:\Google Drive\MATLAB\PhD\ColorFERET_data_ALL');
    %mat_filename = [pathstr '\..\Data\CUFSF_ColorFERET_allData_VIP_plusNormPhotosAffine_photosSubset_842_old2.mat'];
    mat_filename = ['M:\University\From Home 12 Mar15\Transformation and Recognition Program\Data\CUFSF_ColorFERET_allData_VIP_plusNormPhotosAffine_photosSubset_842.mat']; % USE BELOW LINE in final implementation
    %mat_filename = [pathstr '\..\Data\CUFSF_ColorFERET_allData_VIP_plusNormPhotosAffine_photosSubset_842_new.mat'];
    if ~exist(mat_filename,'file')
        warning('Ensure that database paths are correct in ''CUFSF_read_photo_pts'' and ''CUFSF_read_sketch_pts3'' ');
        CUFSF_read_photo_pts % USE CUFSF_read_photo_pts2b in final implementation
        CUFSF_read_sketch_pts3
%         dataAll_sketches = dataAll;
%         clear dataAll;
%         dataAll_photos = data_photos.dataAll;

        data_photos.norm_photos     = face_alignment_MAIN(data_photos, 'CG');
        data_sketches.norm_sketches = face_alignment_MAIN(data_sketches, 'CG');
        
        [data_photos_subset, data_photos_subset_neg] = CUFSF_colorFERET_filter(data_photos);
        
        data_photos_sketches = [];
        data_photos_sketches.photos_norm = data_photos_subset.norm_photos;
        data_photos_sketches.sketches_norm = data_sketches.norm_sketches;
        data_photos_sketches.personData_photos = data_photos_subset.personData;
        data_photos_sketches.personData_sketches = data_sketches.personData;
        data_photos_sketches.subj_id = data_photos_subset.subj_id;
        
        data_photos_sketches.personData_photos_neg_full = data_photos_subset_neg.personData;
        data_photos_sketches.photos_norm_neg_full = data_photos_subset_neg.norm_photos;
        
        data_photos_subset_neg.subj_id_val = str2double(data_photos_subset_neg.subj_id);
        subset_neg_unique_nos = unique(data_photos_subset_neg.subj_id_val);
        for s=1:length(subset_neg_unique_nos)
            r = find((subset_neg_unique_nos(s)==data_photos_subset_neg.subj_id_val));
            n = randi(length(r), 1);
            data_photos_sketches.personData_photos_neg{s, 1} = data_photos_subset_neg.personData{r(n)};            
            data_photos_sketches.photos_norm_neg{s, 1} = data_photos_subset_neg.norm_photos{r(n)};
            data_photos_sketches.subj_id_val_neg{s, 1} = data_photos_subset_neg.subj_id_val(r(n));
            data_photos_sketches.personNo_neg(s, 1) = r(n);
        end
        save(['M:\University\From Home 12 Mar15\Transformation and Recognition Program\Data\CUFSF_ColorFERET_allData_VIP_plusNormPhotosAffine_photosSubset_842_new']);
    else
        load(mat_filename);%'F:\Google Drive\MATLAB\Eigentransformation\ColorFERET_norm_ALL_and_sketches__renamed2_and_using_CUFSF_file_names2_newNamingConvention2.mat'); %load('D:\Google Drive\MATLAB\Eigentransformation\ColorFERET_norm_ALL');
        %data_photos_sketches = data_sketches_subset; % Temporary
    end
else
    error('This database is not supported yet');
end

%load('F:\Google Drive\MATLAB\PhD\ColorFERET_data_ALL');

% Derive the number of images
Nimgs = size(data_photos_sketches.subj_id, 1)

%% 0.3 Set the seed of the global random number generator
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);


%% 0.4 Set parameters/preferences
param.transformation_mode = 'face-to-sketch'; %'sketch-to-face';
param.transformation_fn = 'Eigentransformation';%'none';%'Eigentransformation';
param.FR_fn = 'Eigenfaces';
param.filter_gender = 0;
param.filter_ethnicity = 0;
choice = 0; % Determien if compute quality for images

% Select number of images used for training/testing, set1/set2
param.M = 420; %round(0.5*Nimgs); %600; % Number of images for testing
param.N = Nimgs - param.M;
disp(param.M)
disp(param.N)

if ~strcmpi(param.transformation_mode, 'face-to-sketch') && ~strcmpi(param.transformation_mode, 'sketch-to-face')
    error('Incorrect transformation mode');
end
if ~strcmpi(param.transformation_fn, 'Eigentransformation') && ~strcmpi(param.transformation_fn, 'Eigenpatches') && ~strcmpi(param.transformation_fn, 'none')
    error('Incorrect transformation method');
end
if ~strcmpi(param.FR_fn, 'Eigenfaces')
    error('Incorrect face recogniser');
end
if ((param.filter_gender~=0) && (param.filter_gender~=1))
    error('Parameter for gender filtering must be either ''1'' (on) or ''0'' (off)');
end
if param.filter_ethnicity~=0 && param.filter_ethnicity~=1
    error('Parameter for ethnicity filtering must be either ''1'' (on) or ''0'' (off)');
end

%% 1.1 Split data for training and testing of transformation algorithm

[param.Nrows, param.Ncols] = size(data_photos_sketches.photos_norm{1}); % Extract size based on one photo...assuming all photos/sketches imported already have same dimesions

[imgs_train_synth, imgs_train_FR, imgs_test, param.nos, param.imgs_subj_id] = split_photos_sketches_v2(data_photos_sketches, param);

%% 1.2 Convert images from 2-D to 1-D
%param.Nrows = param.Nrows - 15;
%param.Ncols = param.Ncols - 15;
% Pre-allocation
imgs_train_synth.photos_reshaped = zeros(param.Nrows*param.Ncols, length(param.nos.nos_train_synth));
imgs_train_FR.photos_reshaped = zeros(param.Nrows*param.Ncols, length(param.nos.nos_train_FR));
imgs_test.photos_reshaped = zeros(param.Nrows*param.Ncols, length(param.nos.nos_test));
imgs_train_synth.sketches_reshaped = zeros(param.Nrows*param.Ncols, length(param.nos.nos_train_synth));
imgs_train_FR.sketches_reshaped = zeros(param.Nrows*param.Ncols, length(param.nos.nos_train_FR));
imgs_test.sketches_reshaped = zeros(param.Nrows*param.Ncols, length(param.nos.nos_test));

warning('Ensure that reshaping an H x W picture into (H*W) x 1'); % H = Height; W = Width

imgs_train_synth_gender = cell(1, length(param.nos.nos_train_synth));
imgs_train_synth_ethnicity = cell(1, length(param.nos.nos_train_synth));
for n = 1:length(param.nos.nos_train_synth)
    imgs_train_synth.photos_reshaped(:, n) = imgs_train_synth.photos{n}(:); %reshape(imgs1_norm_train{n}, param.Nrows*param.Ncols, 1);
    imgs_train_synth.sketches_reshaped(:, n) = imgs_train_synth.sketches{n}(:); 
    
    imgs_train_synth_gender{n} = data_photos_sketches.personData_photos{param.nos.nos_train_synth(n)}.Gender;
    imgs_train_synth_ethnicity{n} = data_photos_sketches.personData_photos{param.nos.nos_train_synth(n)}.Race;
end

imgs_train_FR_gender = cell(1, length(param.nos.nos_train_FR));
imgs_train_FR_ethnicity = cell(1, length(param.nos.nos_train_FR));
for n = 1:length(param.nos.nos_train_FR)
    imgs_train_FR.photos_reshaped(:, n) = imgs_train_FR.photos{n}(:); %reshape(imgs1_norm_train{n}, param.Nrows*param.Ncols, 1);
    imgs_train_FR.sketches_reshaped(:, n) = imgs_train_FR.sketches{n}(:); 
    
    imgs_train_FR_gender{n} = data_photos_sketches.personData_photos{param.nos.nos_train_FR(n)}.Gender;
    imgs_train_FR_ethnicity{n} = data_photos_sketches.personData_photos{param.nos.nos_train_FR(n)}.Race;
end

imgs_test_gender = cell(1, length(param.nos.nos_test));
imgs_test_ethnicity = cell(1, length(param.nos.nos_test));
for n = 1:length(param.nos.nos_test)
    imgs_test.photos_reshaped(:, n) = imgs_test.photos{n}(:); %reshape(imgs1_norm_train{n}, param.Nrows*param.Ncols, 1);
    imgs_test.sketches_reshaped(:, n) = imgs_test.sketches{n}(:); 
    
    imgs_test_gender{n} = data_photos_sketches.personData_photos{param.nos.nos_test(n)}.Gender;
    imgs_test_ethnicity{n} = data_photos_sketches.personData_photos{param.nos.nos_test(n)}.Race;
end

if strcmpi(param.transformation_mode, 'face-to-sketch')
    imgs1_train_synth = imgs_train_synth.photos_reshaped;
    imgs2_train_synth = imgs_train_synth.sketches_reshaped;
    
    imgs1_train_FR = imgs_train_FR.photos_reshaped;
    imgs2_train_FR = imgs_train_FR.sketches_reshaped;
    
    imgs1_test = imgs_test.photos_reshaped;
    imgs2_test = imgs_test.sketches_reshaped;
    
elseif strcmpi(param.transformation_mode, 'sketch-to-face')
    imgs2_train_synth = imgs_train_synth.photos_reshaped;
    imgs1_train_synth = imgs_train_synth.sketches_reshaped;
    
    imgs2_train_FR = imgs_train_FR.photos_reshaped;
    imgs1_train_FR = imgs_train_FR.sketches_reshaped;
    
    imgs2_test = imgs_test.photos_reshaped;
    imgs1_test = imgs_test.sketches_reshaped;
else
    error('This type of transformation is not supported');
end


%%  REMOVE this section: 1.3 Select gallery, FR training set for sketch-to-face; for face-to-sketch, gallery is synthesized sketches (i.e. photos transformed to sketches)
% 
% % Gallery
% gallery = data_photos_sketches.photos_norm_neg;
% gallery_reshaped = imgs_norm_neg_reshaped;
% gallery_subj_id = cell2mat(data_photos_sketches.subj_id_val_neg);  
% gallery_gender = imgs_neg_gender;
% gallery_ethnicity = imgs_neg_ethnicity;
% 
% % % Select data to use for FR training:
% % if strcmpi(param.transformation_mode, 'sketch-to-face')
% %     imgs_FR_train = gallery;    
% %     imgs_FR_train_reshaped = gallery_reshaped;
% %     imgs_FR_train_subj_id = gallery_subj_id;    
% %     imgs_FR_train_gender = gallery_gender;
% %     imgs_FR_train_ethnicity = gallery_ethnicity;    
% % end
% 
% 
% % Therefore there are 10 matrices (5 types) which will be used:
% % imgs1_norm_train AND imgs1_norm_train_reshaped
% % imgs2_norm_train AND imgs2_norm_train_reshaped
% % imgs1_norm_test AND imgs1_norm_test_reshaped
% % imgs2_norm_test AND imgs2_norm_test_reshaped
% % imgs_FR_train AND imgs_FR_train_reshaped
% 
% 


%% 2.1 Do training for transformation

if strcmpi(param.transformation_fn, 'Eigentransformation')
    % Eigentransformation:
    train_data.imgs_train = imgs1_train_synth; %imgs1_norm_train_reshaped; % Training based on set 1...i.e. will be transforming from image space of set 1 to that of set 2
    train_data.Nimgs = length(param.nos.nos_train_synth);

    mdl_Transformation_train = transformation_fn(train_data, param.transformation_fn, 'training');

elseif strcmpi(param.transformation_fn, 'Eigenpatches')
    % Eigenpatches
    train_data.patch_size_lr = 128;
    train_data.patch_size_hr = train_data.patch_size_lr;
    train_data.overlap_lr = train_data.patch_size_lr/2;
    train_data.overlap_hr = train_data.overlap_lr;
    if strcmpi(param.transformation_mode, 'face-to-sketch')
        train_data.imgs1_train = imgs_train_synth.photos;
        train_data.imgs1_train_reshaped = imgs_train_synth.photos_reshaped;
        train_data.imgs2_train = imgs_train_synth.sketches;
    elseif strcmpi(param.transformation_mode, 'sketch-to-face')
        train_data.imgs1_train = imgs_train_synth.sketches;
        train_data.imgs1_train_reshaped = imgs_train_synth.sketches_reshaped;
        train_data.imgs2_train = imgs_train_synth.photos;
    end
    train_data.Nimgs = length(param.nos.nos_train_synth);
    mdl_Transformation_train = eigentransformation_patch(train_data, 'training');
end

%% 2.2 Set parameters for transformation (i.e. create data set and set plotFlag to determine if plot transformed photos or not

transformation_data = mdl_Transformation_train;

transformation_data.imgs2_train_reshaped = imgs2_train_synth;% imgs2_norm_train_reshaped; %norm_imgs_sketches_train{pic}; % Used for transformation

transformation_data.plotFlag = 0;
if transformation_data.plotFlag==1
    transformation_data.imgs1_train_reshaped = imgs1_train_synth; %imgs1_norm_train_reshaped;
    transformation_data.Nrows = param.Nrows;
    transformation_data.Ncols = param.Ncols;
end

% For eigenpatches:
transformation_data.Nrows = param.Nrows;
transformation_data.Ncols = param.Ncols;

%% 2.3. Transformation        % OLD (ignore) -> Do testing - transform images and do identification

% Compare 'transformed_imgs' (probe) with images in gallery (i.e. [imgs2_train_reshaped imgs2_test_reshaped])


%[res, c] = transformation_rec_test(test_data, train_data_FR, transformation_data, param, imgs1_norm_test_reshaped, imgs1_norm_test, imgs2_norm_test, imgs_test_subj_id);




% Set parameters for image quality evaluation:
eigentransformation_distance = zeros(param.N, 2);
%FSIM = zeros(param.N, 1);
%FSIMc = zeros(param.N, 1);
% MS-SSIM parameters
K = [0.01 0.03];
winsize = 11;
sigma = 1.5;
window = fspecial('gaussian', winsize, sigma);
level = 5;
weight = [0.0448 0.2856 0.3001 0.2363 0.1333];
method = 'product';




%Need to transform pictures in 2 sets: 
% A: imgs_train_FR 
%and 
% B: imgs_test

% A. Transform pictures which will be used to train FR
% Pictures to transform:
if strcmpi(param.transformation_mode, 'sketch-to-face')
    pics_to_transform = imgs_train_FR.sketches; %imgs1_norm_test;
    
elseif strcmpi(param.transformation_mode, 'face-to-sketch')
    pics_to_transform = imgs_train_FR.photos; %gallery;
    
end
pics_to_transform_reshaped  = imgs1_train_FR; %imgs1_norm_test_reshaped;
pics_to_transform_subj_id   = param.imgs_subj_id.imgs_train_FR_subj_id; %imgs_test_subj_id;
pics_to_transform_gender    = imgs_train_FR_gender; %imgs_test_gender;
pics_to_transform_ethnicity = imgs_train_FR_ethnicity; %imgs_test_ethnicity; 

    
transformed_imgs = cell(length(pics_to_transform), 1); %

timeAll.tElapsed_transform = zeros(1, param.N);

transformed_imgs_reshaped = zeros(param.Nrows*param.Ncols, length(pics_to_transform));
for pic = 1:length(pics_to_transform) % % N % For now, use training images for test...so that all images in testing set are in training set and therefore everyone should be recognised (if have enough images for training...otherwise eigenfaces and reconstruction will be inaccurate)
    
    if strcmpi(param.transformation_fn, 'Eigentransformation');    
        
        % CHECK!:
        if transformation_data.plotFlag==1
            transformation_data.img1 = pics_to_transform{pic}; %probe{pic}==imgs1_norm_test{pic}; %
            pics_to_transformNo = pics_to_transform_subj_id(pic);
            transformation_data.img2 = imgs_FR_train{find(imgs_FR_train_subj_id==pics_to_transformNo)}; %imgs2_norm_test{pic}; % If choice==2, img2 will be wrong %
        end
        
        % Set picture to be transformed
        transformation_data.inImg = pics_to_transform_reshaped(:,pic); % % To be transformed
        
        % Get transformed picture
        tStart_transform = tic;
        transformed_img = transformation_fn(transformation_data, param.transformation_fn, 'transformation');
        timeAll.tElapsed_transform(pic) = toc(tStart_transform);
        imgs1_train_FR_transformed{pic} = ((reshape(transformed_img, param.Nrows, param.Ncols)));        
        imgs1_train_FR_transformed_reshaped(:, pic) = imgs1_train_FR_transformed{pic}(:); % OR transformed_img ? % OR transformed_imgs_reshaped(:, pic) = transformed_img;
        
%        if choice==1
%             [quality.FSIM(pic), quality.FSIMc(pic)] = FeatureSIM(imgs2_norm_test{pic}, (reshape(transformed_imgs{pic}, param.Nrows, param.Ncols)), 'gray');
%             quality.MSSSIM(pic) = ssim_mscale_new(double(imgs2_norm_test{pic}), double(reshape(transformed_imgs{pic}, param.Nrows, param.Ncols)), K, window, level, weight, 'product');
%             quality.PSNR(pic) = (10*log10((255^2)./mean2((double(imgs2_norm_test{pic}) - double(reshape(transformed_imgs{pic}, param.Nrows, param.Ncols))).^2)));
%        end


    elseif strcmpi(param.transformation_fn, 'Eigenpatches');    
        
        % CHECK!:
        if transformation_data.plotFlag==1
            transformation_data.img1 = pics_to_transform{pic}; %probe{pic}==imgs1_norm_test{pic}; %
            pics_to_transformNo = pics_to_transform_subj_id(pic);
            transformation_data.img2 = imgs_FR_train{find(imgs_FR_train_subj_id==pics_to_transformNo)}; %imgs2_norm_test{pic}; % If choice==2, img2 will be wrong %
        end
        
        % Set picture to be transformed
        transformation_data.inImg = pics_to_transform{pic};%(:,pic); % % To be transformed
        
        % Get transformed picture
        tStart_transform = tic;
        transformed_img = eigentransformation_patch(transformation_data, 'transformation'); %transformation_fn(transformation_data, param.transformation_fn, 'transformation');
        timeAll.tElapsed_transform(pic) = toc(tStart_transform);
        imgs1_train_FR_transformed{pic} = ((reshape(transformed_img, param.Nrows, param.Ncols)));        
        imgs1_train_FR_transformed_reshaped(:, pic) = imgs1_train_FR_transformed{pic}(:); % OR transformed_img ? % OR transformed_imgs_reshaped(:, pic) = transformed_img;
        



        

    % CHECK! :
    elseif strcmpi(param.transformation_fn, 'None')
        FR_test_data.inImg_reshaped = pics_to_transform_reshaped(:,pic); %transformed_imgs{pic}; %reshape(transformed_imgs{pic}, Nrows, Ncols); %norm_imgs_train{pic};
    end
    %figure; subplot(1,2,1); imshow(mat2gray(pics_to_transform{pic})); subplot(1,2,2); imshow(mat2gray(imgs1_train_FR_transformed{pic}));
end

% B. Transform pictures which will be used for testing
% Pictures to transform:
if strcmpi(param.transformation_mode, 'sketch-to-face')
    pics_to_transform = imgs_test.sketches; %imgs1_norm_test; 
    orig_transformed_pics = imgs_test.photos;
elseif strcmpi(param.transformation_mode, 'face-to-sketch')
    pics_to_transform = imgs_test.photos; %gallery;
    orig_transformed_pics = imgs_test.sketches;
end
pics_to_transform_reshaped  = imgs1_test; %imgs1_norm_test_reshaped;
pics_to_transform_subj_id   = param.imgs_subj_id.imgs_test_subj_id; %imgs_test_subj_id;
pics_to_transform_gender    = imgs_test_gender; %imgs_test_gender;
pics_to_transform_ethnicity = imgs_test_ethnicity; %imgs_test_ethnicity; 

    
transformed_imgs = cell(length(pics_to_transform), 1); %

eigentransformation_distance = zeros(param.N, 2);

timeAll.tElapsed_transform = zeros(1, param.N);

transformed_imgs_reshaped = zeros(param.Nrows*param.Ncols, length(pics_to_transform));

n = 0;
fprintf('\n');
for pic = 1:length(pics_to_transform) % % N % For now, use training images for test...so that all images in testing set are in training set and therefore everyone should be recognised (if have enough images for training...otherwise eigenfaces and reconstruction will be inaccurate)
    
    msg = sprintf('Picture %d/%d of testing set', pic, length(pics_to_transform));
    fprintf(repmat('\b', 1, n));
    fprintf(msg);
    n = numel(msg);

    if strcmpi(param.transformation_fn, 'Eigentransformation');    
        
        % CHECK!:
        if transformation_data.plotFlag==1
            transformation_data.img1 = pics_to_transform{pic}; %probe{pic}==imgs1_norm_test{pic}; %
            pics_to_transformNo = pics_to_transform_subj_id(pic);
            transformation_data.img2 = imgs_FR_train{find(imgs_FR_train_subj_id==pics_to_transformNo)}; %imgs2_norm_test{pic}; % If choice==2, img2 will be wrong %
        end
        
        % Set picture to be transformed
        transformation_data.inImg = pics_to_transform_reshaped(:,pic); % % To be transformed
        
        % Get transformed picture
        tStart_transform = tic;
        transformed_img = transformation_fn(transformation_data, param.transformation_fn, 'transformation');
        timeAll.tElapsed_transform(pic) = toc(tStart_transform);
        transformed_imgs{pic} = ((reshape(transformed_img, param.Nrows, param.Ncols)));        
        transformed_imgs_reshaped(:, pic) = transformed_imgs{pic}(:); % OR transformed_img ? % OR transformed_imgs_reshaped(:, pic) = transformed_img;
        
%        if choice==1
%             [quality.FSIM(pic), quality.FSIMc(pic)] = FeatureSIM(imgs2_norm_test{pic}, (reshape(transformed_imgs{pic}, param.Nrows, param.Ncols)), 'gray');
%             quality.MSSSIM(pic) = ssim_mscale_new(double(imgs2_norm_test{pic}), double(reshape(transformed_imgs{pic}, param.Nrows, param.Ncols)), K, window, level, weight, 'product');
%             quality.PSNR(pic) = (10*log10((255^2)./mean2((double(imgs2_norm_test{pic}) - double(reshape(transformed_imgs{pic}, param.Nrows, param.Ncols))).^2)));
%        end
       if choice==1
            [quality.FSIM(pic), quality.FSIMc(pic)] = FeatureSIM(orig_transformed_pics{pic}, transformed_imgs{pic}, 'gray');
            quality.MSSSIM(pic) = ssim_mscale_new(double(orig_transformed_pics{pic}), double(transformed_imgs{pic}), K, window, level, weight, 'product');
            quality.PSNR(pic) = (10*log10((255^2)./mean2((double(orig_transformed_pics{pic}) - double(reshape(transformed_imgs{pic}, param.Nrows, param.Ncols))).^2))); %((10*log10((255^2)))./(mean2((double(orig_transformed_pics{pic}) - double(transformed_imgs{pic}))).^2));
       end
       
    elseif strcmpi(param.transformation_fn, 'Eigenpatches');
        
        % CHECK!:
        if transformation_data.plotFlag==1
            transformation_data.img1 = pics_to_transform{pic}; %probe{pic}==imgs1_norm_test{pic}; %
            pics_to_transformNo = pics_to_transform_subj_id(pic);
            transformation_data.img2 = imgs_FR_train{find(imgs_FR_train_subj_id==pics_to_transformNo)}; %imgs2_norm_test{pic}; % If choice==2, img2 will be wrong %
        end
        
        % Set picture to be transformed
        transformation_data.inImg = pics_to_transform{pic};%(:,pic); % % To be transformed
        
        % Get transformed picture
        tStart_transform = tic;
        transformed_img = eigentransformation_patch(transformation_data, 'transformation'); %transformation_fn(transformation_data, param.transformation_fn, 'transformation');
        timeAll.tElapsed_transform(pic) = toc(tStart_transform);
        transformed_imgs{pic} = ((reshape(transformed_img, param.Nrows, param.Ncols)));
        transformed_imgs_reshaped(:, pic) = transformed_imgs{pic}(:); % OR transformed_img ? % OR transformed_imgs_reshaped(:, pic) = transformed_img;
        
        if choice==1
            [quality.FSIM(pic), quality.FSIMc(pic)] = FeatureSIM(orig_transformed_pics{pic}, transformed_imgs{pic}, 'gray');
            quality.MSSSIM(pic) = ssim_mscale_new(double(orig_transformed_pics{pic}), double(transformed_imgs{pic}), K, window, level, weight, 'product');
            quality.PSNR(pic) = (10*log10((255^2)./mean2((double(orig_transformed_pics{pic}) - double(reshape(transformed_imgs{pic}, param.Nrows, param.Ncols))).^2))); %((10*log10((255^2)))./(mean2((double(orig_transformed_pics{pic}) - double(transformed_imgs{pic}))).^2));
        end
        
        % CHECK! :
    elseif strcmpi(param.transformation_fn, 'None')
        %FR_test_data.inImg_reshaped = pics_to_transform_reshaped(:,pic); %transformed_imgs{pic}; %reshape(transformed_imgs{pic}, Nrows, Ncols); %norm_imgs_train{pic};
        transformed_imgs{pic} = pics_to_transform{pic};
        transformed_imgs_reshaped(:, pic) = pics_to_transform_reshaped(:,pic);
    end
    %figure; subplot(1,3,1); imshow(mat2gray(pics_to_transform{pic})); subplot(1,3,2); imshow(mat2gray(orig_transformed_pics{pic})); subplot(1,3,3); imshow(mat2gray(transformed_imgs{pic}));
end

% % % % % % % % % % % % % % % % % % % %

%% 3. Do training for face identification      % Old [ignore] --> ...using all photos/sketches in database
% E.g. Normal eigenfaces

% 3.1 Select data to use for FR training:
if strcmpi(param.transformation_fn, 'none') 
    imgs_FR_train =  [imgs_train_FR.photos]; %transformed_imgs;    
    imgs_FR_train_reshaped = [imgs_train_FR.photos_reshaped]; %transformed_imgs_reshaped; 
elseif strcmpi(param.transformation_mode, 'sketch-to-face')
    imgs_FR_train =  [imgs_train_FR.photos; imgs1_train_FR_transformed']; %transformed_imgs;    
    imgs_FR_train_reshaped = [imgs_train_FR.photos_reshaped imgs1_train_FR_transformed_reshaped]; %transformed_imgs_reshaped; 
elseif strcmpi(param.transformation_mode, 'face-to-sketch')
    imgs_FR_train =  [imgs_train_FR.sketches; imgs1_train_FR_transformed']; %transformed_imgs;    
    imgs_FR_train_reshaped = [imgs_train_FR.sketches_reshaped imgs1_train_FR_transformed_reshaped]; %transformed_imgs_reshaped; 
end
imgs_FR_train_subj_id = [param.imgs_subj_id.imgs_train_FR_subj_id; param.imgs_subj_id.imgs_train_FR_subj_id]; %gallery_subj_id;    
imgs_FR_train_gender = [imgs_train_FR_gender imgs_train_FR_gender]; %gallery_gender;
imgs_FR_train_ethnicity = [imgs_train_FR_ethnicity imgs_train_FR_ethnicity]; %gallery_ethnicity; 

% 3.2 Do training
train_data_FR.imgs_train = imgs_FR_train_reshaped; %[imgs2_test_reshaped]; % [imgs2_train_reshaped imgs2_test_reshaped];    %imgs2_norm_train;
train_data_FR.nos = (imgs_FR_train_subj_id); %param.nos_test;
train_data_FR.Nimgs = size(train_data_FR.imgs_train, 2); %M+N; % M+N=Nimgs

mdl_FR_train = Eigenfaces5(train_data_FR, 'training'); %face_rec(train_data_FR, param.FR_fn, 'training'); % For eigenfaces: If photo-to-sketch, have sketch eigenspace; if sketch-to-photo, have photo eigenspace

%% 3.3 TEST
% Set parameters for FR
% FR_test_data = mdl_FR_train;
% FR_test_data.M = param.M;
% FR_test_data.nos = (imgs_FR_train_subj_id); %param.nos_test;


% FR test:

% PROBE, GALLERY:
if strcmpi(param.transformation_mode, 'sketch-to-face')
    probe = transformed_imgs;
    probe_reshaped = transformed_imgs_reshaped;
    
    gallery = imgs_test.photos;
    gallery_reshaped = imgs_test.photos_reshaped;  
    
    bp = zeros(length(mdl_FR_train.Evalue), length(gallery));
    br = zeros(length(mdl_FR_train.Evalue), length(probe));    
elseif strcmpi(param.transformation_mode, 'face-to-sketch')
    probe = imgs_test.sketches; %imgs2_norm_test;
    probe_reshaped = imgs_test.sketches_reshaped; %imgs2_norm_test_reshaped;
    
    gallery = transformed_imgs;
    gallery_reshaped = transformed_imgs_reshaped;
    
    br = zeros(length(mdl_FR_train.Evalue), length(gallery));
    bs = zeros(length(mdl_FR_train.Evalue), length(probe));    
end
% Probe and Gallery have same subj-ID/Gender/Ethnicity (since both come from test set):
probe_subj_id = param.imgs_subj_id.imgs_test_subj_id;
probe_gender = imgs_test_gender;
probe_ethnicity = imgs_test_ethnicity;

gallery_subj_id = param.imgs_subj_id.imgs_test_subj_id;
gallery_gender = imgs_test_gender;
gallery_ethnicity = imgs_test_ethnicity;

res = cell(1, length(probe));
c = zeros(length(probe), 1);

% Weight projection of probes and gallery AND comparison
if strcmpi(param.transformation_mode, 'sketch-to-face')
    et_s2p = load('F:\Google Drive\University\PhD\Logs\042 - 12Mar15_Home_newEigenfacesCont\Matlab_12Mar15_T48__withQuality2', 'quality_FSIMc', 'quality', 'quality_PSNR', 'quality_MSSSIM', 'Nranks', 'res_rank', 'param', 'res');
    %HAOG = load('F:\Google Drive\University\PhD\Logs\042 - 12Mar15_Home_newEigenfacesCont\Matlab_12Mar15_T57', 'Nranks', 'res_rank', 'param', 'res');
    HAOG = load('F:\Google Drive\University\PhD\Logs\053 - 28Mar15 - Fusion cont\Matlab_28Mar15_T57b.mat', 'Nranks', 'res_rank', 'param', 'res');

    for pic = 1:length(gallery) % length(param.nos.nos_test)
        mdl_FR_train.inImg_reshaped = gallery_reshaped(:, pic);
        bp(:, pic) = Eigenfaces5(mdl_FR_train, 'testing2'); % Gallery weights
    end
    for pic = 1:length(probe) % length(param.nos.nos_test)
        mdl_FR_train.inImg_reshaped = probe_reshaped(:, pic);
        br(:, pic) = Eigenfaces5(mdl_FR_train, 'testing2'); % Probe weights (of pseudo-photos)
        picNo_actual = param.nos.nos_test(pic);
        
        theta = bp - repmat(br(:, pic), [1,size(bp,2)]);        
        %res{pic} = sqrt(sum(theta.^2, 1)); %theta4 
        
        p=2;
        res{pic} = (sum(abs(theta).^p, 1)).^(1/p);    
%         res{pic} = (res{pic}-min(res{pic}))/(max(res{pic})-min(res{pic}));
%         et_s2p.res_orig{pic}(:, et_s2p.res{pic}(2, :)) = et_s2p.res{pic};
%         et_s2p.res_orig{pic}(1, :) = (et_s2p.res_orig{pic}(1, :)-min(et_s2p.res_orig{pic}(1, :)))/(max(et_s2p.res_orig{pic}(1, :))-min(et_s2p.res_orig{pic}(1, :)));
%         res{pic} = res{pic} + et_s2p.res_orig{pic}(1, :);
%         HAOG.res_orig{pic}(:, HAOG.res{pic}(2, :)) = HAOG.res{pic};
%         HAOG.res_orig{pic}(1, :) = (HAOG.res_orig{pic}(1, :)-min(HAOG.res_orig{pic}(1, :)))/(max(HAOG.res_orig{pic}(1, :))-min(HAOG.res_orig{pic}(1, :)));
%         res{pic} = res{pic} + HAOG.res_orig{pic}(1, :);
        
%         for i=1:length(param.nos.nos_test)
%             res{pic}(1, i) = norm(theta(:,i),2);%sqrt(sum(theta.^2, 1)); %theta4
%         end
        
        % Filtering using gender information:
        if param.filter_gender==1
            picNo_actual_gender = probe_gender(pic);
            gender_filter = find(~strcmpi(picNo_actual_gender, gallery_gender));
            res{pic}(1, gender_filter) = inf;
            filter.gender_sum(pic) = length(gender_filter);
        end
        % Filtering using ethnicity information:
        if param.filter_ethnicity==1
            picNo_actual_ethnicity = probe_ethnicity(pic);
            ethnicity_filter = find(~strcmpi(picNo_actual_ethnicity, gallery_ethnicity));
            res{pic}(1, ethnicity_filter) = inf;
            filter.ethnicity_sum(pic) = length(ethnicity_filter);
        end        
        res{pic}(3, :) = param.nos.nos_test'; %FR_test_data.nos;
        [res{pic}(1, :), res{pic}(2, :)] = sort(res{pic}(1, :), 'ascend');
        res{pic}(3,:) = res{pic}(3, res{pic}(2, :));
        c(pic) = find(res{pic}(3,:)==picNo_actual); % [~,c(pic),~] = 
    end  
elseif strcmpi(param.transformation_mode, 'face-to-sketch')
%    et_p2s = load('F:\Google Drive\University\PhD\Logs\042 - 12Mar15_Home_newEigenfacesCont\Matlab_12Mar15_T49__withQuality3', 'quality_FSIMc', 'quality', 'quality_PSNR', 'quality_MSSSIM', 'Nranks', 'res_rank', 'param', 'res');
    %HAOG = load('F:\Google Drive\University\PhD\Logs\042 - 12Mar15_Home_newEigenfacesCont\Matlab_12Mar15_T57', 'Nranks', 'res_rank', 'param', 'res');
%    HAOG = load('F:\Google Drive\University\PhD\Logs\053 - 28Mar15 - Fusion cont\Matlab_28Mar15_T57b.mat', 'Nranks', 'res_rank', 'param', 'res');

    for pic = 1:length(gallery)
        mdl_FR_train.inImg_reshaped = gallery_reshaped(:, pic);
        br(:, pic) = Eigenfaces5(mdl_FR_train, 'testing2'); % Gallery weights (of pseudo-sketches)
    end
    for pic = 1:length(probe)
        mdl_FR_train.inImg_reshaped = probe_reshaped(:, pic);
        bs(:, pic) = Eigenfaces5(mdl_FR_train, 'testing2'); % Probe weights
        picNo_actual = param.nos.nos_test(pic);
        
        theta = br - repmat(bs(:, pic), [1,size(br,2)]);
        %res{pic} = sqrt(sum(theta.^2, 1)); %theta4        
        
        p=2;
        res{pic} = (sum(abs(theta).^p, 1)).^(1/p);    
%         res{pic} = (res{pic}-min(res{pic}))/(max(res{pic})-min(res{pic}));
%         et_p2s.res_orig{pic}(:, et_p2s.res{pic}(2, :)) = et_p2s.res{pic};
%         et_p2s.res_orig{pic}(1, :) = (et_p2s.res_orig{pic}(1, :)-min(et_p2s.res_orig{pic}(1, :)))/(max(et_p2s.res_orig{pic}(1, :))-min(et_p2s.res_orig{pic}(1, :)));
%         res{pic} = res{pic} + et_p2s.res_orig{pic}(1, :);
%         HAOG.res_orig{pic}(:, HAOG.res{pic}(2, :)) = HAOG.res{pic};
%         HAOG.res_orig{pic}(1, :) = (HAOG.res_orig{pic}(1, :)-min(HAOG.res_orig{pic}(1, :)))/(max(HAOG.res_orig{pic}(1, :))-min(HAOG.res_orig{pic}(1, :)));
%         res{pic} = res{pic} + HAOG.res_orig{pic}(1, :);
%
%        res{pic} = (1.0.*((sum(abs(theta).^p, 1)).^(1/p)))+(1.0.*(et_p2s.res_orig{pic}(1, :)));
%
%         for i=1:length(param.nos.nos_test)
%             res{pic}(1, i) = norm(theta(:,i),2);%sqrt(sum(theta.^2, 1)); %theta4
%         end

        
        % Filtering using gender information:
        if param.filter_gender==1
            picNo_actual_gender = probe_gender(pic);
            gender_filter = find(~strcmpi(picNo_actual_gender, gallery_gender));
            res{pic}(1, gender_filter) = inf;
            filter.gender_sum(pic) = length(gender_filter);
        end
        % Filtering using ethnicity information:
        if param.filter_ethnicity==1
            picNo_actual_ethnicity = probe_ethnicity(pic);
            ethnicity_filter = find(~strcmpi(picNo_actual_ethnicity, gallery_ethnicity));
            res{pic}(1, ethnicity_filter) = inf;
            filter.ethnicity_sum(pic) = length(ethnicity_filter);
        end
        res{pic}(3, :) = param.nos.nos_test'; %FR_test_data.nos;
        [res{pic}(1, :), res{pic}(2, :)] = sort(res{pic}(1, :), 'ascend');
        res{pic}(3,:) = res{pic}(3, res{pic}(2, :));
        c(pic) = find(res{pic}(3,:)==picNo_actual); % [~,c(pic),~] = 
    end
end


% % % % % % % % % % % % % % % % % % % %



% res = cell(1, param.N);
% 
% c = zeros(param.N, 1);
% timeAll.tElapsed_FR = zeros(1, param.N);
% for pic = 1:length(probe) %param.N % % N % For now, use training images for test...so that all images in testing set are in training set and therefore everyone should be recognised (if have enough images for training...otherwise eigenfaces and reconstruction will be inaccurate)
%     
%     FR_test_data.inImg_reshaped = probe_reshaped(:, pic); % Current probe image to compare with gallery %reshape(transformed_imgs{pic}, Nrows, Ncols); %norm_imgs_train{pic};
%         
%     % Do recognition
%     tStart_FR = tic;
%     res{pic} = face_rec(FR_test_data, param.FR_fn, 'testing'); % Sorted results
%     timeAll.tElapsed_FR(pic) = toc(tStart_FR);
%     if strcmpi(param.FR_fn,'Eigenfaces')
%         picNo_actual = probe_subj_id(pic);
%         
%         % Filtering using gender information:
%         if param.filter_gender==1
%             picNo_actual_gender = probe_gender(pic);
%             gender_filter = find(~strcmpi(picNo_actual_gender, gallery_gender));
%             res{pic}(1, gender_filter) = inf;
%             filter.gender_sum(pic) = length(gender_filter);
%         end
%         % Filtering using ethnicity information:
%         if param.filter_ethnicity==1
%             picNo_actual_ethnicity = probe_ethnicity(pic);
%             ethnicity_filter = find(~strcmpi(picNo_actual_ethnicity, gallery_ethnicity));
%             res{pic}(1, ethnicity_filter) = inf;
%             filter.ethnicity_sum(pic) = length(ethnicity_filter);
%         end
%         
%         res{pic}(3, :) = FR_test_data.nos;
%         [res{pic}(1, :), res{pic}(2, :)] = sort(res{pic}(1, :), 'ascend');
%         res{pic}(3,:) = res{pic}(3, res{pic}(2, :));
%         
%   %      [res{pic}(1, :), res{pic}(2, :)] = sort(res{pic}, 'ascend');    
%   %      res{pic}(3,:) = train_data_FR.nos(res{pic}(2, :));  % Actual picture numbers of sorted gallery images 
%         
%         %res{pic}(4, :) = imgs_test_subj_id(res{pic}(2, :));
%         %picNo_actual = probe_subj_id(pic); %imgs_test_subj_id(pic);   %train_data_FR.nos(pic); %==nos_test(pic) % Actual picture number of test image
%         c(pic) = find(res{pic}(3,:)==picNo_actual); % [~,c(pic),~] = 
%         eigentransformation_distance(pic, 1) = res{pic}(1, c(pic));
%         eigentransformation_distance(pic, 2) = res{pic}(1, 1);
%     end
% end
% eigentransformation_distance(end+2, :) = mean(eigentransformation_distance, 1);
% 
% timeAll.tElapsed_transform(end+2) = mean(timeAll.tElapsed_transform);
% timeAll.tElapsed_FR(end+2) = mean(timeAll.tElapsed_FR);
% if param.filter_gender==1
%     filter.gender_sum(end+2) = sum(filter.gender_sum);
%     disp(filter.gender_sum(end))
% end
% if param.filter_ethnicity==1
%     filter.ethnicity_sum(end+2) = sum(filter.ethnicity_sum);
%     disp(filter.ethnicity_sum(end))
% end

fprintf('\n');

% 5. Get Results
% Rank-n identification accuracy
Nranks = length(param.nos.nos_test); % length(imgs_FR_train_subj_id); %param.N; % 188
res_rank = zeros(Nranks, 2);

for rank = 1:Nranks
    res_rank(rank, 1) = sum(c<=rank); % Number of images
    res_rank(rank, 2) = (res_rank(rank, 1)./length(param.nos.nos_test))*100; % Percentage of all images
    if rank<=10 % Display first 10 ranks
        fprintf('Rank %03d: %0.4f\n', rank, res_rank(rank, 2));
    end
end

if strcmpi(param.transformation_fn, 'Eigentransformation')
    figure('Name', 'rank-n accuracy'); plot(1:Nranks, res_rank(:, 2), '-x'); xlabel('Rank'); ylabel('Percentage of correctly classified faces'); title(['Identification vs. rank for M=' num2str(param.M) ' and ' num2str(size(mdl_Transformation_train.w,1)) ' eigenvectors used'])
else
    figure('Name', 'rank-n accuracy'); plot(1:Nranks, res_rank(:, 2), '-x'); xlabel('Rank'); ylabel('Percentage of correctly classified faces'); title(['Identification vs. rank for M=' num2str(param.M)])
end
grid on

if choice==1 && ~strcmpi(param.transformation_fn, 'none');
    % Quality using FSIM/FSIMc
    ctr = 1;
    for q = 0:0.01:1
        quality_FSIMc(ctr, 1) = sum(quality.FSIMc<=q); % Number of images
        quality_FSIMc(ctr, 2) = (quality_FSIMc(ctr, 1)./length(param.nos.nos_test))*100; % Percentage of all images
        
        quality_MSSSIM(ctr, 1) = sum(quality.MSSSIM<=q); % Number of images
        quality_MSSSIM(ctr, 2) = (quality_MSSSIM(ctr, 1)./length(param.nos.nos_test))*100; % Percentage of all images
        
%         quality_PSNR(ctr, 1) = sum(quality.PSNR<=q); % Number of images
%         quality_PSNR(ctr, 2) = (quality_PSNR(ctr, 1)./length(param.nos.nos_test))*100; % Percentage of all images
        ctr = ctr + 1;
    end
    figure('Name', 'Quality between transformed images and corresponding original image using FSIMc and MS-SSIM'); plot(0:0.01:1, quality_FSIMc(:, 2), '-x', 0:0.01:1, quality_MSSSIM(:, 2), '-rx'); xlabel('Quality'); ylabel('Percentage of pictures having measured quality'); title(['FSIM for M=' num2str(param.M) ' and N=' num2str(param.N)]); legend('FSIMc', 'MS-SSIM');
    grid on;

    rho_p_FSIMc = corr(c, quality.FSIMc')
    rho_s_FSIMc = corr(c, quality.FSIMc', 'type', 'Spearman')
    
    rho_p_MSSSIM = corr(c, quality.MSSSIM')
    rho_s_MSSSIM = corr(c, quality.MSSSIM', 'type', 'Spearman')
    
    ctr = 1;
    for q = 0:0.01:ceil(max(quality.PSNR))        
        quality_PSNR(ctr, 1) = sum(quality.PSNR<=q); % Number of images
        quality_PSNR(ctr, 2) = (quality_PSNR(ctr, 1)./length(param.nos.nos_test))*100; % Percentage of all images
        ctr = ctr + 1;
    end
    figure('Name', 'Quality between transformed images and corresponding original image using PSNR'); plot(0:0.01:ceil(max(quality.PSNR)), quality_PSNR(:, 2), '-x');
    grid on;

    rho_p_PSNR = corr(c, quality.PSNR')
    rho_s_PSNR = corr(c, quality.PSNR', 'type', 'Spearman')

    fprintf('Mean FSIMc = %0.4f +/- %0.4f standard deviations\n', mean(quality.FSIMc), std(quality.FSIMc));
    fprintf('Mean MS-SSIM = %0.4f +/- %0.4f standard deviations\n', mean(quality.MSSSIM), std(quality.MSSSIM));
    fprintf('Mean PSNR = %0.4f +/- %0.4f standard deviations\n', mean(quality.PSNR), std(quality.PSNR));
   
end

%% Delete unnecessary variables

%clearvars -except CLK tStart imgs_FR_train_subj_id choice timeAll filter param res c res_rank mdl_Transformation_train mdl_FR_train quality eigentransformation_distance Nranks;

%%
% for ii=1:length(imgs_test_subj_id_sorted)
%     if sum(imgs_test_subj_id_sorted(ii)==imgs_FR_train_subj_id)==0
%         fprintf('ii=%d',ii);
%     end
% end
%
% T006 vs T007 vs T014 vs T018:
% T06 = load('F:\Google Drive\PhD\Logs\013 - 29Dec14 - Tests\Matlab_29Dec14_T006', 'res_rank');
% T14 = load('F:\Google Drive\PhD\Logs\013 - 29Dec14 - Tests\Matlab_29Dec14_T014', 'res_rank');
% T18 = load('F:\Google Drive\PhD\Logs\013 - 29Dec14 - Tests\Matlab_29Dec14_T018', 'res_rank');
% T07 = load('F:\Google Drive\PhD\Logs\013 - 29Dec14 - Tests\Matlab_29Dec14_T007', 'res_rank');
% figure; plot(1:168, T06.res_rank(:, 2), 1:168, T07.res_rank(:, 2),  1:168, T14.res_rank(:, 2),   1:168, T18.res_rank(:, 2)); legend('Face-to-sketch [old alignment]', 'Face-to-sketch [new alignment]', 'Sketch-to-face [new alignment]', 'No transformation [new alignment] [Gallery=photos, Probe=sketches]'); grid on;

%% Subplots to show test image, corresponding image in other domain and best matches
% 
% test_img = 4;
% test_img_no = param.nos_test(test_img);
% 
% Ngallery = length(imgs_FR_train);
% Nranks_plot = 5;
% 
% 
% 
% figure('Name', ['Probe/gallery/transformation comparison for test_img ' num2str(test_img)]);
% if strcmpi(param.transformation_mode, 'face-to-sketch') % VIP: imgs_FR_train uses imgs2_norm_test (sketches)...corresponding photos are in img1_norm_test
%     warning('Ensure imgs1_all is correct (i.e. depends on images used for imgs_FR_train (gallery))');
%     imgs1_all = [imgs1_norm_test; imgs1_norm_train]; % Depends on what images used for 'imgs_FR_train'
% 
%     subplot(2,Nranks_plot+3, 1); imshow(probe{test_img}); title('Probe image');
%     subplot(2,Nranks_plot+3, 2); imshow(mat2gray(reshape(transformed_imgs{test_img}, param.Nrows,param.Ncols))); title('Tramnsformed Probe image');
%     subplot(2,Nranks_plot+3, 3); imshow(gallery{res{test_img}(2,c(test_img))}); title('Matching photo in gallery');
%     for rank=1:Nranks_plot
%         subplot(2,Nranks_plot+3,3+rank); imshow(gallery{res{test_img}(2,rank)}); title(['Rank-' num2str(rank)]);
%         subplot(2,Nranks_plot+3, Nranks_plot+3+3+rank); imshow(imgs1_all{res{test_img}(2,rank)}); title(['Rank-' num2str(rank)]);
%     end
% elseif strcmpi(param.transformation_mode, 'sketch-to-face')
%     subplot(1,Nranks_plot+3, 1); imshow(probe{test_img}); title('Probe image');
%     subplot(1,Nranks_plot+3, 2); imshow(mat2gray(reshape(transformed_imgs{test_img}, param.Nrows,param.Ncols))); title('Tramnsformed Probe image');
%     subplot(1,Nranks_plot+3, 3); imshow(gallery{res{test_img}(2,c(test_img))}); title('Matching photo in gallery');
%     for rank=1:Nranks_plot
%         subplot(1,Nranks_plot+3,3+rank); imshow(gallery{res{test_img}(2,rank)}); title(['Rank-' num2str(rank)]);
%     end
% elseif strcmpi(param.transformation_mode, 'none')
% else
%     error('Incorrect transformation mode');
% end



%% END
fprintf('\nReady!\n\n');
% save rred_live_recap rred 

X = load('gong');
player = audioplayer(X.y, X.Fs);
play(player)

tElapsed = toc(tStart); % End timer
CLK2 = clock;
fprintf('Function ended at %d/%d/%d, %d:%02.0f:%02.0f \n\n', CLK2(3),CLK2(2),CLK2(1),CLK2(4),CLK2(5),CLK2(6) );
 
fprintf('Function completed in %2.2f seconds\n', tElapsed);
if ( CLK2(4) < CLK(4) ) % For midnight
    hours = (24+CLK2(4))-CLK(4);
elseif ( (CLK2(4) < CLK(4)) && (CLK2(5)<CLK(5)) )
    hours = (24+CLK2(4))-CLK(4) -1;
elseif ( (CLK2(4) > CLK(4)) && (CLK2(5)<CLK(5)) ) % Change of hour dows not mean 1 hour
    hours = CLK2(4)-CLK(4) -1;
else
    hours = CLK2(4)-CLK(4);
end
if ((CLK2(5) < CLK(5)) && CLK2(6) < CLK(6)) % change in minutes but no. of seconds of end time smaller than that of start time ex. 19:10:43, 19:18:33
    mins = 60+CLK2(5) - CLK(5) -1;
elseif CLK2(5) < CLK(5)
    mins = (60+CLK2(5))-CLK(5);
elseif ((CLK2(5) > CLK(5)) && CLK2(6) < CLK(6))
    mins = CLK2(5) - CLK(5) -1;
else
    mins = CLK2(5)-CLK(5);
end
if CLK2(6) < CLK(6)
    secs = (60+CLK2(6))-CLK(6);
else
    secs = CLK2(6)-CLK(6);
end
fprintf('i.e. Function duration: %d:%02.0f:%02.0f \n\n', hours, mins, secs );

%send_mail_message('chris_p_galea@yahoo.co.uk', 'MATLAB ready', [mfilename ' is ready']);
 
%profile viewer;
 
% if(tElapsed > 30)
%     error('myApp:Check','WARNING: \nFunction has taken too long!') %Low-quality image!\n
% end
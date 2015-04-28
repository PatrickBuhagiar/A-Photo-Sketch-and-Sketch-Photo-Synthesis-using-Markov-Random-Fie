function [Training_Sketch_Patches, Training_Photo_Patches] = InitializeTrainingSet_search_Patch_training(train_sketch_path, train_photo_path, sketch_list, image_list, model_path, n_images, patch_size, overlap_size, Input_Patches)

%% Initialisation
load([model_path 'AOM_MultiPIE_InTheWild']);
Training_Sketches = cell(1,n_images);
Training_Photos = Training_Sketches;
Training_Sketch_Patches = cell(1,n_images);
Training_Photo_Patches = Training_Sketch_Patches;
%% Iterate Images

for i = 1:n_images,
    %% Read Skecth
    filename = strcat(train_sketch_path, sketch_list(i).name);
    Train_In = imread(filename);
   
    % Define input image
    [~, ~, n_channels] = size(Train_In); 
    if n_channels == 3
      Train_In = double(rgb2gray(Train_In))./255;
    else
      Train_In = double(Train_In)./255;
    end 

    % Run Matlab's standard face detector
    faceDetector = vision.CascadeObjectDetector();
    bbox = step(faceDetector, Train_In);
        
    %Initialize parameters
    level = 1;
    interpolation = {'none', 'bilinear'};
    param.AAM.interpolation = cell2mat(interpolation(2));
    over_scale = 0.9 * bbox(end,3) / AAM.resolution{1}(1);
    current_shape = AAM.shape_mean_level{level} * over_scale - repmat(mean(AAM.shape_mean_level{level} * over_scale, 1), AAM.n_vertices, 1) +  repmat([(bbox(end,1)+bbox(end,3)/2), ((bbox(end,2)+bbox(end,4)/1.65))], AAM.n_vertices, 1) / 2^(level-1);
    current_shape2 = current_shape * 2^(level-1);
    
    %Warp image
    [warped_img] = warpImage(AAM.shape_mean_scaled{level}, AAM.texture_base{level}, AAM.triangles, AAM.resolution{level}, current_shape2, Train_In, param.AAM.interpolation);
    Training_Sketches{i} = SubSample(warped_img,67, 67);
    %figure;imshow(Training_Sketches{i});
    
    %% Read photo
    filename = strcat(train_photo_path, image_list(i).name);
    Train_Out = imread(filename);
   
    % Define input image
    [~, ~, n_channels] = size(Train_Out); 
    if n_channels == 3
      Train_Out = double(rgb2gray(Train_Out))./255;
    else
      Train_Out = double(Train_Out)./255;
    end 
    
    % Run Matlab's standard face detector
    faceDetector = vision.CascadeObjectDetector();
    bbox = step(faceDetector, Train_Out);
        
    %Initialize parameters
    level = 1;
    interpolation = {'none', 'bilinear'};
    param.AAM.interpolation = cell2mat(interpolation(2));
    over_scale = 0.9 * bbox(end,3) / AAM.resolution{1}(1);
    current_shape = AAM.shape_mean_level{level} * over_scale - repmat(mean(AAM.shape_mean_level{level} * over_scale, 1), AAM.n_vertices, 1) +  repmat([(bbox(end,1)+bbox(end,3)/2), ((bbox(end,2)+bbox(end,4)/1.65))], AAM.n_vertices, 1) / 2^(level-1);
    current_shape2 = current_shape * 2^(level-1);
    
    %Warp image
    [warped_img] = warpImage(AAM.shape_mean_scaled{level}, AAM.texture_base{level}, AAM.triangles, AAM.resolution{level}, current_shape2, Train_Out, param.AAM.interpolation);
    Training_Photos{i} = SubSample(warped_img,67, 67);
    %figure;imshow(Training_Photos{i});
    
    
    %% create patches (DEFINE PATCH SIZE AND OVERLAP SIZE HERE)
    
    
    [Training_Sketch_patches, Training_Photo_patches] = patches_search_Patch(patch_size, overlap_size, Training_Sketches{i}, Training_Photos{i}, Input_Patches{i});
    Training_Sketch_Patches{i} = Training_Sketch_patches;
    Training_Photo_Patches{i} = Training_Photo_patches;

end

end

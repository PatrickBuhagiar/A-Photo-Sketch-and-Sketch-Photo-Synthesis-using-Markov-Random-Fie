%clc; close all; clear all;

sigma = 1;
patch_size = 10;
overlap_size = 5;
K_candidates = 5;
full_image_size = [67 67];

%% acquire TrainingSketches path

train_sketch_path = '..\02_data\CUHK_Sketches\';
sketch_list = dir([train_sketch_path '*jpg']);
n_sketches = size(sketch_list, 1);


%% acquire TrainingImages path

train_photo_path = '..\02_data\CUHK_Photos\';
image_list = dir([train_photo_path '*jpg']);
n_images = size(image_list, 1);

%% acquire Input Image path
test_path = '..\02_data\03_InTheWild\';
input_list = dir([test_path '*jpg']);
n_input = size(input_list, 1);

%% Load Model for warp fit
addpath('\01_functions\');
model_path = '..\02_data\02_Multi-PIE\01_AOMs\';
% load([model_path 'AOM_MultiPIE_InTheWild']);

%% Obtain Patches - TO CHANGE
%Each one contains a cell array with 3D arrays. This 3D array consists of 
%all the patches and are indexed by the third dimension. the length of the 
%cell array is equal to the number of images.  

%Training_Sketch_Patches = InitializeTrainingSet(train_sketch_path, sketch_list, model_path, n_sketches, patch_size, overlap_size);
%Training_Photo_Patches = InitializeTrainingSet(train_photo_path, image_list, model_path, n_images,  patch_size, overlap_size);
Input_Patches = InitializeTrainingSet_search_Patch(test_path, input_list, model_path, n_input,  patch_size, overlap_size);
[Training_Sketch_Patches, Training_Photo_Patches] = InitializeTrainingSet_search_Patch_training(train_sketch_path, train_photo_path, sketch_list, image_list, model_path, n_images, patch_size, overlap_size, Input_Patches); 
n_patches = length(Input_Patches{1});



%% Remains the same

% %% Obtain best matches
candidate_patches = PseudoImage(Input_Patches, Training_Sketch_Patches, Training_Photo_Patches, K_candidates);

% calculate phi
phi_patches = Calculate_phi(Input_Patches, candidate_patches, sigma, K_candidates);

%calculate neighbours
neighbours = Determine_Neighbours(n_patches, patch_size, overlap_size, full_image_size);

%decide patches
final_patches_index = decide_patches(candidate_patches, phi_patches, neighbours, overlap_size, full_image_size);


%extract chosen patches
final_patches = zeros(patch_size, patch_size, n_patches);
for i=1:n_patches,
    final_patches(:,:,i) = candidate_patches(:,:,final_patches_index(1,i),i);
end

I = GeneratePsuedoPhoto(final_patches, full_image_size, overlap_size);
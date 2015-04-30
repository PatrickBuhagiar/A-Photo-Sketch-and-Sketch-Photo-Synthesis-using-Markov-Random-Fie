%%CRITERIA

%DATASET: 88 testing, 94 training

%PATCH SIZE: OVERLAPS
%    5: 4,3,2
%   10: 5,3,2 
%   15: 7,5,3
%   20: 10,5,3

%CANDIDATES: 5,10,15

%Quality determined by SSIM and PSNR
%RESULT: SSIM and PSNR scores for 94 images.

%First sketch to photo, then vice versa

%% Parameters
sigma = 1;
patch_sizes = [20 15 10 5;10 7 5 2];
K_candidates = [5];
full_image_size = [160 160];

%Load Model for warp fit
addpath('\01_functions\');
model_path = '..\02_data\02_Multi-PIE\01_AOMs\';

%% Load Dataset

%acquire Training Sketches path
train_sketch_path = '..\02_data\New_CUHK_Training_Sketches\';
sketch_list = dir([train_sketch_path '*jpg']);
n_training_sketches = size(sketch_list, 1);

%acquire Training Photos path
train_photo_path = '..\02_data\New_CUHK_Training_Photos\';
image_list = dir([train_photo_path '*jpg']);
n__training_photos = size(image_list, 1);

%acquire Testing Sketches path
test_sketch_path = '..\02_data\New_CUHK_Testing_Sketches\';
input__skecth_list = dir([test_sketch_path '*jpg']);
n_sketch_input = size(input__skecth_list, 1);

%acquire Testing Photos path
test_photo_path = '..\02_data\New_CUHK_Testing_Photos\';
input__photo_list = dir([test_photo_path '*jpg']);
n_photo_input = size(input__photo_list, 1);

%% Loops
[n_ol,n_ps] = size(patch_sizes);
 mkdir('Search_Path_Optimisation');
 
for ps = 1:n_ps,
    patch_size = patch_sizes(1,ps);
    for ol = 2: n_ol,
        overlap_size = patch_sizes(ol, ps);
        for kc = 1: length(K_candidates),
            K =  K_candidates(kc);
            patch_size
            overlap_size
            K
            %% acquire patches
            [Testing_Sketch_Patches, Testing_Sketch_Warped] = InitializeTrainingSet(test_sketch_path, input__skecth_list, model_path, n_sketch_input,  patch_size, overlap_size);
            [Testing_Photo_Patches, Testing_Photo_Warped] = InitializeTrainingSet(test_photo_path, input__photo_list, model_path, n_photo_input,  patch_size, overlap_size);
             
            %The difference lies here. choosing training patches depend on testing patches
            [Training_Sketch_Patches, Training_Photo_Patches] = InitializeTrainingSet_search_Patch_training(train_sketch_path, train_photo_path, sketch_list, image_list, model_path, n__training_photos, patch_size, overlap_size, Testing_Sketch_Patches);
            
            
            
            n_patches = length(Testing_Photo_Patches{1});
            fileID = fopen('Search_Path_Optimisation_test.txt','a+');
            fprintf(fileID,'\nfor %d PATCH SIZE, %d OVERLAP SIZE, %d CANDIDATES:\n',patch_size, overlap_size, K);
            fclose(fileID);
            
            ssim_A = zeros(1, n_photo_input);
            ssim_B = zeros(1, n_photo_input);
            psnr_A = zeros(1, n_photo_input);
            psnr_B = zeros(1, n_photo_input);
            for t_i=1:n_photo_input, %for each training image
                
                %% Obtain best matches
                disp('Obtaining best matches.');
                candidate_patches_A = PseudoImage(Testing_Sketch_Patches{t_i}, Training_Sketch_Patches, Training_Photo_Patches, K);
                candidate_patches_B = PseudoImage(Testing_Photo_Patches{t_i}, Training_Photo_Patches, Training_Sketch_Patches, K);
                
                %% calculate phi
                disp('Calculating Phi');
                phi_patches_A = Calculate_phi(Testing_Sketch_Patches{t_i}, candidate_patches_A, sigma, K);
                phi_patches_B = Calculate_phi(Testing_Photo_Patches{t_i}, candidate_patches_B, sigma, K);
                
                %% calculate neighbours
                disp('Calculating Neighbours');
                neighbours = Determine_Neighbours(n_patches, patch_size, overlap_size, full_image_size);
                
                %% decide patches
                 disp('Deciding Patches');
                final_patches_index_A = decide_patches(candidate_patches_A, phi_patches_A, neighbours, overlap_size, full_image_size);
                final_patches_index_B = decide_patches(candidate_patches_B, phi_patches_B, neighbours, overlap_size, full_image_size);
                
                %% extract chosen patches
                disp('Extracting Chosen Patches');
                final_patches_A = zeros(patch_size, patch_size, n_patches);
                final_patches_B = zeros(patch_size, patch_size, n_patches);
                for i=1:n_patches,
                    final_patches_A(:,:,i) = candidate_patches_A(:,:,final_patches_index_A(1,i),i);
                    final_patches_B(:,:,i) = candidate_patches_B(:,:,final_patches_index_B(1,i),i);
                end
                
                %% Generate Image
                disp('Generating Images');
                I_A = GeneratePsuedoPhoto(final_patches_A, full_image_size, overlap_size);
                I_B = GeneratePsuedoPhoto(final_patches_B, full_image_size, overlap_size);
                %figure; imshow(I_A); figure; imshow(I_B);
                mkdir(['images_A_' num2str(patch_size) '_' num2str(overlap_size) '_' num2str(K) ]);
                mkdir(['images_B_' num2str(patch_size) '_' num2str(overlap_size) '_' num2str(K) ]);
                I_A_file = ['images_A_' num2str(patch_size) '_' num2str(overlap_size) '_' num2str(K) '\I_A_'  num2str(t_i) '.jpg'];
                I_B_file = ['images_B_' num2str(patch_size) '_' num2str(overlap_size) '_' num2str(K) '\I_B_'  num2str(t_i) '.jpg'];
                imwrite(I_A, I_A_file);
                imwrite(I_B, I_B_file);
                
                [x,y] = size(Testing_Sketch_Warped{t_i});

                %Resize image to fit patches
                Xstart = 1:(patch_size - overlap_size):x;
                Ystart = 1:(patch_size - overlap_size):y;   
                Xend = round(patch_size:(patch_size - overlap_size):x);
                Yend = round(patch_size:(patch_size - overlap_size):y); 
                Input_A = SubSample( Testing_Sketch_Warped{t_i},Xend(end),Yend(end));
                Input_B = SubSample(Testing_Photo_Warped{t_i},Xend(end),Yend(end));

                PSNR_A = psnr(I_A, Input_A);
                PSNR_B = psnr(I_B, Input_B);
                
                SSIM_A = ssim(I_A, Input_A);
                SSIM_B = ssim(I_B, Input_B);
                
                psnr_A(t_i) = PSNR_A;
                psnr_B(t_i) = PSNR_B;
                ssim_A(t_i) = SSIM_A;
                ssim_B(t_i) = SSIM_B;
                
                %% Export to file
                disp('Exporting to file');
                fileID = fopen('Sketch_Photo_Test_A.txt', 'a+');
                fprintf(fileID,'For test images %d, PSNR_A: %.2f | PSN_B: %.2f | SSIM_A: %.2f | SSIM_B: %.2f \n',t_i, PSNR_A, PSNR_B, SSIM_A, SSIM_B);
                fclose(fileID);
            end
            %% csv write
            disp('CSV write results');
            ssim_A_file = ['Search_Path_Optimisation\ssimcomp_ssim_A_' num2str(patch_size) '_' num2str(overlap_size) '_' num2str(K) 'I_A_'  num2str(t_i) '.csv'];
            ssim_B_file = ['Search_Path_Optimisation\ssimcomp_ssim_B_' num2str(patch_size) '_' num2str(overlap_size) '_' num2str(K) 'I_B_'  num2str(t_i) '.csv'];
            psnr_A_file = ['Search_Path_Optimisation\ssimcomp_psnr_A_' num2str(patch_size) '_' num2str(overlap_size) '_' num2str(K) 'I_A_'  num2str(t_i) '.csv'];
            psnr_B_file = ['Search_Path_Optimisation\ssimcomp_psnr_B_' num2str(patch_size) '_' num2str(overlap_size) '_' num2str(K) 'I_B_'  num2str(t_i) '.csv'];
            csvwrite(ssim_A_file, ssim_A);
            csvwrite(ssim_B_file, ssim_B);
            csvwrite(psnr_A_file, psnr_A);
            csvwrite(psnr_B_file, psnr_B);
        end     
    end   
end









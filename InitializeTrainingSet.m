function Patches = InitializeTrainingSet(train_sketch_path, sketch_list, model_path, n_sketches, patch_size, overlap_size),

%% Initialisation
load([model_path 'AOM_MultiPIE_InTheWild']);
Patches = cell(1,n_sketches);

%% Iterate Images

for i = 1:n_sketches,
   
    %Read Image 
    filename = strcat(train_sketch_path, sketch_list(i).name);
    I = imread(filename);
   
    % Define input image
    [~, ~, n_channels] = size(I); 
    if n_channels == 3
      I = double(rgb2gray(I))./255;
    else
      I = double(I)./255;
    end 

    % Run Matlab's standard face detector
    faceDetector = vision.CascadeObjectDetector();
    bbox = step(faceDetector, I);
        
    %Initialize parameters
    level = 1;
    interpolation = {'none', 'bilinear'};
    param.AAM.interpolation = cell2mat(interpolation(2));
    over_scale = 0.9 * bbox(end,3) / AAM.resolution{1}(1);
    current_shape = AAM.shape_mean_level{level} * over_scale - repmat(mean(AAM.shape_mean_level{level} * over_scale, 1), AAM.n_vertices, 1) +  repmat([(bbox(end,1)+bbox(end,3)/2), ((bbox(end,2)+bbox(end,4)/1.65))], AAM.n_vertices, 1) / 2^(level-1);
    current_shape2 = current_shape * 2^(level-1);
    
    %Warp image
    [warped_img] = warpImage(AAM.shape_mean_scaled{level}, AAM.texture_base{level}, AAM.triangles, AAM.resolution{level}, current_shape2, I, param.AAM.interpolation);
    %figure;imshow(warped_img);
    %create patches (DEFINE PATCH SIZE AND OVERLAP SIZE HERE)
    Patches{1,i} = patches(warped_img, patch_size, overlap_size);
   
end

end

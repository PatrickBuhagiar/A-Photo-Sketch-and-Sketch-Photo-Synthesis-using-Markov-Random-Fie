train_path = '..\02_data\New_CUHK_Training_Photos\';
train_list =  dir([train_path '*jpg']);
n_images = size(train_list, 1);

model_path = '..\02_data\02_Multi-PIE\01_AOMs\';
load([model_path 'AOM_MultiPIE_InTheWild']);

mkdir('photos_warped_images3');
mkdir('photos_real_images3');

for j=1:n_images
       j
       filename = strcat(train_path, train_list(j).name);
       X = double(imread(filename))./255;
       % Run Matlab's standard face detector
       faceDetector = vision.CascadeObjectDetector();
       bbox = step(faceDetector, X);

       %Initialize parameters
       level = 1;
       interpolation = {'none', 'bilinear'};
       param.AAM.interpolation = cell2mat(interpolation(2));
       over_scale = 0.9 * bbox(end,3) / AAM.resolution{1}(1);
       current_shape = AAM.shape_mean_level{level} * over_scale - repmat(mean(AAM.shape_mean_level{level} * over_scale, 1), AAM.n_vertices, 1) +  repmat([(bbox(end,1)+bbox(end,3)/2), ((bbox(end,2)+bbox(end,4)/1.65))], AAM.n_vertices, 1) / 2^(level-1);
       current_shape2 = current_shape * 2^(level-1);

       %Warp image
       warped_img = warpImage(AAM.shape_mean_scaled{level}, AAM.texture_base{level}, AAM.triangles, AAM.resolution{level}, current_shape2, X, param.AAM.interpolation);
       write = ['photos_warped_images3\image_'  num2str(100+j) '.jpg'];
       imwrite(warped_img, write);
       write = ['photos_real_images3\image_'  num2str(100+j) '.jpg'];
       imwrite(X, write);

end


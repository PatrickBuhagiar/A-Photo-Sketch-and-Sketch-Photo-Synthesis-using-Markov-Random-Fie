close all
clear all
clc
cd(fileparts(which(mfilename)));
cd ..

path = '.\02_data\01_XM2VTS\';


%% Read Images

test_path = [path '02_test_images\'];
load([test_path 'image_list.mat']);
load([test_path 'shape_list.mat']);
n_data = size(image_list, 1);


%% Load Model

model_path = [path '\01_AOMs\'];
load([model_path 'AOM_XM2VTS']);


%% Specify Control Parameters

verbose = {true, false};
param.verbose = cell2mat(verbose(1));
param.AAM.verbose = param.verbose;

display = {true, false};
param.display = cell2mat(display(1));
param.AAM.display = param.display;

save_image = {true, false};
param.save_image = cell2mat(save_image(2));

video = {true, false};
param.video = cell2mat(video(1));


%% Fit

ground = zeros(68, 2, n_data);
fitted = zeros(68, 2, n_data);
for i = 1:n_data 
  
  try
    image_path = [test_path image_list(i).name];
    shape_path = [test_path shape_list(i).name];
    ground(:,:,i) = readShape(shape_path, 68);
    disp(['Fitting image: ' int2str(i) '/' int2str(n_data)]);
    [fitted(:,:,i)] = POICA_AOM_XM2VTS(param, AAM, image_path); 
  catch fitting_err
    disp('Fitting_Error!');
  end
    
end

plotFinalResults4(ground, fitted);


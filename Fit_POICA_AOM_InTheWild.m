close all
clear all
clc

cd(fileparts(which(mfilename)));
cd ..

addpath('01_code\01_functions\');

%% Read Images

test_path = '.\02_data\03_InTheWild\';
image_list = dir([test_path '*png']);
n_images = size(image_list, 1);


%% Load Model

model_path = '.\02_data\02_Multi-PIE\01_AOMs\';
load([model_path 'AOM_MultiPIE_InTheWild']);


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
param.video = cell2mat(video(2));


%% Fit

fitted = zeros(68, 2, n_images);
for i = 1:n_images 
  
  try
    image_path = [test_path image_list(i).name];
    disp(['Fitting image: ' int2str(i) '/' int2str(n_images)]);
    [fitted(:,:,i)] = POICA_AOM_InTheWild(param, AAM, image_path); 
  catch fitting_err
    disp('Fitting_Error!');
  end
    
end


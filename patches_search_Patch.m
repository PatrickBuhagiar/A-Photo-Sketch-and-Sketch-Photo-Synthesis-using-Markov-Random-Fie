function [Training_Sketch_Patches, Training_Photo_Patches] = patches_search_Patch(patch_size, overlap_size, Train_Sketch, Train_Photo, Input_Patches)
%patch size 
%overlap size
%Train_in = 67x67
%Train_out = 67x67 (Train_in) pair
%input_patches = XxYxn_patches
n_patches = size(Input_Patches, 3); 
Training_Sketch_Patches = zeros(patch_size, patch_size, n_patches);
Training_Photo_Patches = Training_Sketch_Patches;

[x,y] = size(Train_Sketch);

%Resize image to fit patches
Xstart = 1:(patch_size - overlap_size):x;
Ystart = 1:(patch_size - overlap_size):y;   
Xend = round(patch_size:(patch_size - overlap_size):x);
Yend = round(patch_size:(patch_size - overlap_size):y);

%search space indices
X_s_start = Xstart(1:size(Xend,2)) - 2;
X_s_start(1) = 1;

X_s_end = Xend(1:size(Xend,2)) + 2;
X_s_end(end) = Xend(end);
%temp = SubSample(imIn,Xend(end),Yend(end));

height = length(Xend);
width = length(Yend);


numberOfPatches = size(Xend,2) * size(Yend,2);
patches = zeros(patch_size, patch_size, numberOfPatches); 
c= 1;

for a=1:size(Xend,2),
    for b=1:size(Yend,2),
        test_patch = Input_Patches(:,:,c);
        sketch_search_space = Train_Sketch(X_s_start(a):1: X_s_end(a), X_s_start(b):1: X_s_end(b));
        photo_search_space = Train_Photo(X_s_start(a):1: X_s_end(a), X_s_start(b):1: X_s_end(b));

        f = im2col(sketch_search_space, [patch_size patch_size], 'sliding');
        f2 = im2col(photo_search_space, [patch_size patch_size], 'sliding');

        H = repmat(test_patch(:), 1, size(f,2));

        g =  sum(abs(f - H));
        [im,index] = min(g);

        Training_Sketch_Patches(:,:,c) = reshape(f(:,index), [patch_size,patch_size]);
        Training_Photo_Patches(:,:,c) = reshape(f2(:,index), [patch_size,patch_size]);
        c = c+1;
    end
end 
end
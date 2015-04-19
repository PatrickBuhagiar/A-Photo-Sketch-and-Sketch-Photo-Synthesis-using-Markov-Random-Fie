function [Training_Sketch_Patches, Training_Photo_Patches] = patches_search_Patch(patchSize, overlapSize, Train_In, Train_out, input_patches)
%patch size 
%overlap size
%Train_in = 67x67
%Train_out = 67x67 (Train_in) pair
%input_patches = XxYxn_patches




[x,y] = size(Train_In);
 
 %Resize image to fit patches
 Xstart = 1:(patchSize - overlapSize):x;
 Ystart = 1:(patchSize - overlapSize):y;   
 Xend = round(patchSize:(patchSize - overlapSize):x);
 Yend = round(patchSize:(patchSize - overlapSize):y); 
 %temp = SubSample(imIn,Xend(end),Yend(end));
 
 height = length(Xend);
 width = length(Yend);
 
 
 numberOfPatches = size(Xend,2) * size(Yend,2);
 patches = zeros(patchSize, patchSize, numberOfPatches); 
 c= 1;
 
 for a=1:size(Xend,2),
    for b=1:size(Yend,2),
        patches(:,:,c) = temp(Xstart(a):1:Xend(a),Ystart(b):1:Yend(b),:);
        c= c+ 1;    
    end
 end 
end
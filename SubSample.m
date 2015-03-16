function imOut = SubSample(imIn, height, width)
%Obtain Dimensions 
[x,y,z] = size(imIn);

%initialise output image
imOut = zeros(height, width, z);

%Obtain Ratio
M = x/height;
N = y/width;

%Assignment
imOut(:,:,:) = imIn(round(1:M:x),round(1:N:y), :);

end
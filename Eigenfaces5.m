%function model = eigentransformation_training(Ih_list, xl)

% Based on Eigenfaces1.m
% This is a full implementation of Eigenfaces algorithm
%
% 'Training':
% 'inData' must contain fields for 'M' (number of images used for training), 'norm_imgs_train' (images used for training)
% 'outData' will contain 'E' (eigenfaces), 'c' (weights for reconstruction), 'w' (weights) and 'm' (mean image obtained from training samples)
%
% 'Testing'
% 'inData' must contain data from training, namely 'E' (eigenfaces) 'm' (mean image obtained from training samples) 'w' (weights) and 'inImg' (input image to identify from amongst pictures in training set)
% 'outData' will contain sorted person matches, in ascending order (i.e. first match is most likely match (rank-1) and so on

function outData = Eigenfaces5(inData, mode)


if strcmpi((mode), ('training'))

    % Rename mat-file and get number of images used for training
    
    norm_imgs_train_reshaped = inData.imgs_train;
    
    Nimgs = inData.Nimgs; %size(norm_imgs_train, 1);
    
    % Derive the dimensions of the images
    %[Nrows, Ncols] = size(norm_imgs_train{1});
    

    % Compute the mean vector
    m = mean(norm_imgs_train_reshaped,2);

    % Derive the difference vector L
    L = norm_imgs_train_reshaped - repmat(m,[1,Nimgs]);

    % Compute the covariance matrix
    C = L'*L;

    % Derive the eigen values and eigenvectors of the covariance matrix
    [Evector,Evalue] = eig(C);    
    
    % Derive the eigenvectors and eigenvalues that cover 99% of the variance
    [Evector,Evalue] = eigenvector_selection(Evector,diag(Evalue),0.99);

    % Derive the Eigenfaces
    E = L * Evector * diag(1./sqrt(Evalue));

%     %USE SVD INSTEAD OF EIG (LIKE USED IN PHD_TOOL):
%     % comment out line deriving 'c'    
%     [U2, V2, L2] = svd(C);
%     E = normc(L*U2(:,1:end-1));
%     Evalue = diag(V2); 


    % Derive the projection parameters (weights)
    w = E'*(norm_imgs_train_reshaped - repmat(m,[1,Nimgs]));

    % Derive the reconstruction weights
    c = Evector * diag(1./sqrt(Evalue))* w;
    

    outData.E = E;
    outData.c = c; % Can remove
    outData.w = w;
    outData.m = m;
    %outData.norm_imgs_train_reshaped = norm_imgs_train_reshaped;    
    outData.Evector = Evector;
    outData.Evalue = Evalue;

elseif strcmpi((mode), ('testing')) % Compare 'img1' with 'img2'
    %% Input image

    % Load data
    E = inData.E;
    m = inData.m;    
    w = inData.w;
    
    inImg_reshaped = inData.inImg_reshaped; 

    w_inImg = E'*(double(inImg_reshaped) - m);

    % Identification:
    theta = w - repmat(w_inImg, [1,size(w,2)]);
    theta4 = sqrt(sum(theta.^2, 1));

%     x = w;
%     covar = inv(cov(x(1:end,:)'));
%     y = repmat(w_inImg, [1,size(w,2)]);
%     norm_x = sqrt(x'*covar*x);
%     norm_y = sqrt(y'*covar*y);
%     d = - (x'*covar*y)/(norm_x*norm_y);
    


%     res = [];
%     [res(1, :), res(2, :)] = sort(theta4, 'ascend');
%     
%     res(3, :) = inData.nos(res(2,:)); % Contains actual image numbers
    
    res = theta4; % %
    
    outData = res;
    
elseif strcmpi((mode), ('testing2')) % Compare 'img1' with 'img2'
    %% Input image

    % Load data
    E = inData.E;
    m = inData.m;    
    w = inData.w;
    
    inImg_reshaped = inData.inImg_reshaped; 

    w_inImg = E'*(double(inImg_reshaped) - m);

    outData = w_inImg;
    
elseif strcmpi((mode), ('transformation'))

    % Load data
    E = inData.E;
    m = inData.m;    
    %w = inData.w;
    %M = inData.M;

    %c = inData.c;
    Evector = inData.Evector;
    Evalue = inData.Evalue;
      
    imgs2 = inData.imgs2_train_reshaped; %reshape_imgs(M, Nrows, Ncols, inData.norm_imgs_sketches_train);
    m2 = mean(imgs2, 2);
    
    inImg_reshaped = inData.inImg;    

    w_inImg = E'*(double(inImg_reshaped) - m);
    c_inImg = Evector * diag(1./sqrt(Evalue))* w_inImg; % Reconstruction weights

    %L2 = imgs2 - repmat(m2,[1,size(imgs2,2)]); % This is the 'mathematically correct' version (used in line below instead of imgs2), but can reduce computation time by using imgs2 directly and still get same results
    inImg_reconstructed  = imgs2*c_inImg + m2;
    
    if inData.plotFlag==1
        % Uncomment these for debugging:
        inImg_reconstructed1 = E*w_inImg + m;    
        imgs1 = inData.imgs1_train_reshaped;  
        inImg_reconstructed2 = imgs1*c_inImg + m; % Should be like inImg_reconstructed1
        Nrows = inData.Nrows;
        Ncols = inData.Ncols;
        img1 = inData.img1;
        img2 = inData.img2;
        figure('Name', 'True photo vs. transformed photo, true sketch vs. transformed sketch');
        subplot(1,4,1); imshow(mat2gray(img1)); title('Original img1');
        subplot(1,4,2); imshow(mat2gray(reshape(inImg_reconstructed2(:, 1), Nrows, Ncols))); title('Transformed img1');
        subplot(1,4,3); imshow(mat2gray(img2)); title('Original img2');
        subplot(1,4,4); imshow(mat2gray(reshape(inImg_reconstructed, Nrows, Ncols))); title('Transformed img2');
    end
    
    outData = inImg_reconstructed;
  
    
else
    error('Unkown mode...value must be either ''training'' or ''testing''');
end

end

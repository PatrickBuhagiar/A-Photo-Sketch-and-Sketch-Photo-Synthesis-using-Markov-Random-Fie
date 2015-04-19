function imOut = PseudoImage(imIn, Training_sketches, Training_Photo_Patches, K)

% imIn : Input Image divided into patches 
% imTraining : Training Sketches divided into patches
% imTrainingPair : Photo pairs of Training Sketches
% imOut : image patches

%code segments:
%true for nearest K candidates
candidates = true;
%% Score nearest patch ONLY
%obtain number of patches
if (~candidates)
    TotalImages = length(Training_sketches);
    TotalPatches = size(imIn,3);
    PatchMatches = zeros(1,TotalPatches); %indexes of closest image patches

    %Iterate all Patches of input image
    for i=1:TotalPatches,
        ReferencePatch = imIn(:,:,i);

        PatchScores = zeros(1, TotalImages);
        for j=1:TotalImages,  
            %Use SSIM that returns ordered list.
            PatchScores(j) = ssim(Training_sketches{j}(:,:,i), ReferencePatch); %sum(sum(sqrt((Training_sketches{j}(:,:,i) - ReferencePatch).^2)));
            %aa = ssim(Training_sketches{j}(:,:,i), ReferencePatch);
        end
        [HighestScore, I] = max(PatchScores);
        PatchMatches(i) = I;
    end

    %group nearest patches

    NearestPatches = zeros(size(imIn,1), size(imIn,2), TotalPatches);

    for i=1: TotalPatches,
        j = PatchMatches(i);
        NearestPatches(:,:,i) = Training_Photo_Patches{1,j}(:,:,i);
    end

    imOut = NearestPatches;

%% Score nearest candidate patches only
else
    TotalImages = length(Training_sketches);
    TotalPatches = size(imIn,3);
    %K = 5; %number of candidates
    K_Candidates_sketches_indexes = zeros(K,TotalPatches); %indexes of closest image patches

    %Iterate all Patches of input image
    for i=1:TotalPatches,
        ReferencePatch = imIn(:,:,i);

        PatchScores = zeros(1, TotalImages);
        for j=1:TotalImages,  
            %Use SSIM that returns ordered list.
            
            PatchScores(j) = ssim(Training_sketches{j}(:,:,i), ReferencePatch);% sum(sum(sqrt((Training_sketches{j}(:,:,i) - ReferencePatch).^2)));
            %aa =ssim(Training_sketches{j}(:,:,i), ReferencePatch);
        end
        [scores, indexes] = sort(PatchScores, 'descend');
        temp = indexes(1:K); %take top K
        K_Candidates_sketches_indexes(:,i) = temp(:); %save results
    end

    %group nearest patches
    %K_Candidates_photos(x,y,K,n_patches). x and y are patch sizes
    K_Candidates_photos = zeros(size(imIn,1), size(imIn,2), K, TotalPatches);

    for i=1: TotalPatches,
        k = K_Candidates_sketches_indexes(:,i);
        for j=1: K,
            K_Candidates_photos (:,:,j,i) = Training_Photo_Patches{1,k(j)}(:,:,i); 
        end
        
    end

    imOut = K_Candidates_photos;
    
end
end
function final_patches = decide_patches(candidate_patches, phi_patches, neighbours, overlap, full_image_size)
%Candidate_Patches = (x,y,k_candidate, n_patch)
%phi_patches = (k_candidate, n_patch)
%neighbours = cell array, length n_patches. each cell contains vector of neighbours
%overlap = overlap size

n_patches = size(candidate_patches,4);
n_candidates = size(phi_patches,1);

final_patches = zeros(1, n_patches); %candidate indexes
argmax = zeros(1,n_candidates);


for n=1:n_patches,
    for k=1:n_candidates,
        %obtain overlapping patches for current n_th patch
        neighbour_overlaps = get_overlaps(candidate_patches(:,:,k,n), candidate_patches, neighbours{n}, n, k, overlap, full_image_size);
        %first patch is current overlap. 
        %for this candidate, find the argmax
        argmax(1,k) = neighbour_argmax(neighbour_overlaps, phi_patches, n, k, neighbours);    
    end
    [M, I] = max(argmax);
    final_patches(1,n) = I; %I is candidate index
    
    
end
end
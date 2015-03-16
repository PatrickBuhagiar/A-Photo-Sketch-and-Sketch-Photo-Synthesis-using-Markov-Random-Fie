function overlaps = get_overlaps(candidate_patch, candidate_patches, neighbours, current_patch, current_candidate, overlap, full_image_size)
%candidate_patch = (x,y)
%candidate_patches = (x,y,k_candidates,n_patches)
%neighbours = vector of neighbour indexes
%current_patch = patch index
%current_candidate = candidate index
%overlap = overlap size

n_neighbours = length(neighbours);
n_candidates = size(candidate_patches,3); 
x = size(candidate_patch,1);
overlaps = zeros(x, overlap, n_candidates + 1, n_neighbours); %+1 because includes current patch

%filter out neighbour patches
neighbour_patches = zeros(x,x,n_candidates, n_neighbours);
for n=1:n_neighbours,
    neighbour_patches(:,:,:,n) = candidate_patches(:,:,:,neighbours(n));
    for k=1:n_candidates,
       [d_current, d_neighbour] = get_overlap(current_patch, candidate_patch, neighbours(n), neighbour_patches(:,:,k,n), overlap, full_image_size);   
       overlaps(:,:,1,n) = d_current;
       overlaps(:,:,k+1,n) = d_neighbour;
    end
end

end


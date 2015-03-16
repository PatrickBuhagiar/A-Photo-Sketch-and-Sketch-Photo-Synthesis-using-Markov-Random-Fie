function phi_patches = Calculate_phi(input_patches, candidate_patches, sigma)
%input patches: Cell array (x,y, n_patches)
%candidate_patches: 4D array (x,y,K_candidates, n_patches)
%sigma: standard deviation
%phi_patches: (k_candidate, n_patch)

n_patches = size(candidate_patches,4);
phi_patches = zeros(size(candidate_patches,3),n_patches);
for i=1:n_patches,
    phi_patches(:,i) = phi(input_patches{1}(:,:,i), candidate_patches(:,:,:,i), sigma);     
end
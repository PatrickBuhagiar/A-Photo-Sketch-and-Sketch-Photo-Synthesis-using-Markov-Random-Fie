function argmax = neighbour_argmax(neighbour_overlaps, phi_patches, n_patch, n_candidate, neighbours)
%neighbour_overlaps = (x, overlap, k_candidate+1, n_neighbour)
%phi_patches = (k_candidate, n_patch)
%n_patch = current index index
%n_candidate = current candidate index

n_neighbours = size(neighbour_overlaps,4);
n_candidates = size(phi_patches,1);
sigma = 1;

argmax = 1;

for n=1:n_neighbours,
    Message = zeros(1, n_candidates);
    for k=1:n_candidates,
        Phi = phi_patches(n_candidate, neighbours{n_patch}(n));
        psi = Psi(neighbour_overlaps(:,:,1,n), neighbour_overlaps(:,:,k+1,n), sigma);
        Message(1, n) = Phi*psi;
    end
    [M, I] = max(Message);
    argmax = argmax*M;
end

end
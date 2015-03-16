function p = phi(input_patch, candidate_patches, sigma)
%sigma = standard deviation
%input_patch = 2D matrix patch (x,y)
%candidate_patches = 3D matrix (x,y,K)

y = repmat(input_patch,[1 1 size(candidate_patches, 3)]);
p = sum(sum(exp(-(abs(candidate_patches - y).^2)./(2.*(sigma.^2)))));
p = reshape(p,5,1); %convert 1x1xK to Kx1 vector

end
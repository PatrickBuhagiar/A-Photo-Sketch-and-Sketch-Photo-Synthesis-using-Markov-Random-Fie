function [final_shape, rhos, alphas] = computeWarpUpdate(current_shape, delta, shape_mean, shape_pc, shape_ev, global_transform, n_transformations, shape_dims, triangles, sorted_triangles, p, n_vertices, regularization)
      
  % computeWarpUpdate Summary of this function goes here
  % Detailed explanation goes here

  % Invert Warp on reference shape using ONLY delta_alpha
  if ~isempty(shape_pc)
    delta_alphas = delta(n_transformations + 1:end);
    delta_shape = - shape_pc * delta_alphas;
    delta_shape = reshape(delta_shape, [], 2);
    new_s0_alpha = shape_mean + delta_shape;
  else
    new_s0_alpha = shape_mean;
  end

  % Invert Warp on reference shape using ONLY delta_rho
  delta_rhos = delta(1:n_transformations);
  delta_similarity = - global_transform * delta_rhos;
  delta_similarity = reshape(delta_similarity, [], 2);
  [A t] = computeSimilarityMatricialForm(shape_mean, delta_similarity, triangles, sorted_triangles);
  new_s0_alpha_rho = new_s0_alpha * A + repmat(t, [n_vertices 1]);

  % Compose delta_shape_similarity with current shape
  delta_shape_similarity = new_s0_alpha_rho - shape_mean;
  new_shape = computeWarpComposition(shape_mean, delta_shape_similarity, current_shape, triangles, sorted_triangles);

  % Project new_shape onto silimarity and shape PCA basis
  [final_shape, rhos, alphas] = computeShapeProjection(new_shape, shape_mean, shape_pc, shape_ev, global_transform, shape_dims, triangles, sorted_triangles, p, n_vertices, regularization);

end
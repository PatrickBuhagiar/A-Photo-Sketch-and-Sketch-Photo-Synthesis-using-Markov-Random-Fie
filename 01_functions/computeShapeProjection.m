function [projected_shape, rhos, alphas] = computeShapeProjection(current_shape, shape_mean, shape_pc, shape_ev, global_transform, shape_dims, triangles, sorted_triangles, p, n_vertices, regularization)
      
  % computeShapeProjection Summary of this function goes here
  % Detailed explanation goes here

  % Project new face onto the global similarity transform basis
  rhos = global_transform' * (current_shape(:) - shape_mean(:));

  % Project new face onto shape model once effect of global similarity 
  % is removed transform
  global_similarity = global_transform * rhos;
  global_similarity = reshape(global_similarity, [], 2);
  [A t] = computeSimilarityMatricialForm(shape_mean, global_similarity, triangles, sorted_triangles);
    
  if ~isempty(shape_pc)
    diff_shape = (current_shape - repmat(t, [n_vertices 1])) / A - shape_mean;
    if regularization
      alphas = (eye(shape_dims) + p * diag([ones(1, shape_dims)] ./ shape_ev(1:shape_dims))) \ shape_pc' * diff_shape(:);
    else
      alphas = shape_pc' * diff_shape(:);
    end
%     M_dist = sqrt(sum((alphas' ./ sqrt(shape_ev(1:shape_dims))).^2));
%     alphas = sqrt(shape_dims) * alphas / M_dist;
    pre_final_shape = shape_mean(:) + shape_pc * alphas;
  else
    pre_final_shape = shape_mean(:);
    alphas = [];
  end
  
  % Compute Final shape
  projected_shape = reshape(pre_final_shape, [], 2) * A + repmat(t, [n_vertices 1]);

end
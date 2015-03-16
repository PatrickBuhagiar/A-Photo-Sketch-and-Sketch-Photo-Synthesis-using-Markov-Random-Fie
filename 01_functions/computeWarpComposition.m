function new_shape = computeWarpComposition(mean_shape_scaled, delta_shape, current_shape, triangles, sorted_triangles)
% This function applies piecewise affine image warping according to the
% composed shape of current and the incremental delta shape.

new_shape = zeros(size(mean_shape_scaled));

for i=1:size(mean_shape_scaled,1)
    tx = mean_shape_scaled(i,1) + delta_shape(i,1);
    ty = mean_shape_scaled(i,2) + delta_shape(i,2);
    
    this_triangle = sorted_triangles{i};
    
    if (isempty(this_triangle))
        continue;
    end
    
    v = zeros(length(this_triangle),2);
    
    for j=1:length(this_triangle)
        trij = this_triangle(j);
        v0x = mean_shape_scaled(triangles(trij,1),1);
        v0y = mean_shape_scaled(triangles(trij,1),2);
        v1x = mean_shape_scaled(triangles(trij,2),1);
        v1y = mean_shape_scaled(triangles(trij,2),2);
        v2x = mean_shape_scaled(triangles(trij,3),1);
        v2y = mean_shape_scaled(triangles(trij,3),2);
        
        vt0 = current_shape(triangles(trij,1),1:2);
        vt1 = current_shape(triangles(trij,2),1:2);
        vt2 = current_shape(triangles(trij,3),1:2);
        
        donominator = (v1x - v0x) * (v2y - v0y) - (v1y - v0y) * (v2x - v0x);
        
        a = ((tx - v0x) * (v2y - v0y) - (ty - v0y) * (v2x - v0x)) / donominator;
        b = ((ty - v0y) * (v1x - v0x) - (tx - v0x) * (v1y - v0y)) / donominator;
        
        v(j,:) = vt0 + a * (vt1 - vt0) + b * (vt2 - vt0);
    end
    new_shape(i,:) = mean(v);    
end

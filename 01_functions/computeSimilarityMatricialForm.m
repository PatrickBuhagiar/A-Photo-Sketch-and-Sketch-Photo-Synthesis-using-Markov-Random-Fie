function [A t] = computeSimilarityMatricialForm(mean_shape_scaled, delta_shape, triangles, sorted_triangles)
% This function applies piecewise affine image warping according to the
% composed shape of current and the incremental delta shape.

current_shape = mean_shape_scaled + delta_shape;

% mat = rand(2);
% t2 = rand(1,2);
% current_shape = mean_shape_scaled * mat + repmat(t2, [58 1]);

for i=1:1
  
%     tx = mean_shape_scaled(i,1);
%     ty = mean_shape_scaled(i,2);
%     
%     ttx = current_shape(i,1);
%     tty = current_shape(i,2);
    
    this_triangle = sorted_triangles{i};
    
    if (isempty(this_triangle))
        continue;
    end
    
    for j=1:1
      
        trij = this_triangle(j);
        v0x = mean_shape_scaled(triangles(trij,1),1);
        v0y = mean_shape_scaled(triangles(trij,1),2);
        v1x = mean_shape_scaled(triangles(trij,2),1);
        v1y = mean_shape_scaled(triangles(trij,2),2);
        v2x = mean_shape_scaled(triangles(trij,3),1);
        v2y = mean_shape_scaled(triangles(trij,3),2);
        
        vt0x = current_shape(triangles(trij,1),1);
        vt0y = current_shape(triangles(trij,1),2);
        vt0 = [vt0x vt0y];
        vt1x = current_shape(triangles(trij,2),1);
        vt1y = current_shape(triangles(trij,2),2);
        vt1 = [vt1x vt1y];
        vt2x = current_shape(triangles(trij,3),1);
        vt2y = current_shape(triangles(trij,3),2);
        vt2 = [vt2x vt2y];
        
        
        denominator = (v1x - v0x) * (v2y - v0y) - (v1y - v0y) * (v2x - v0x);
        
        A = v2y - v0y;
        B = v2x - v0x;
        
        C = v1x - v0x;
        D = v1y - v0y;
        
        F = vt1x - vt0x;
        G = vt2x - vt0x;
        
        H = vt1y - vt0y;
        I = vt2y - vt0y;
        
        
        a1 = vt0x + ((-A * v0x + B * v0y) * F + (-C * v0y + D * v0x) * G) / denominator;
        
        a2 = (A * F - D * G) / denominator;
        
        a3 = (C * G - B * F) / denominator;
        
        
        a4 = vt0y + ((-A * v0x + B * v0y) * H + (-C * v0y + D * v0x) * I) / denominator;
        
        a5 = (A * H - D * I) / denominator;
        
        a6 = (C * I - B * H) / denominator;
        
        
        A = [ a2  a5
              a3  a6 ];
           
        t = [ a1  a4 ];
        
    end
    
end

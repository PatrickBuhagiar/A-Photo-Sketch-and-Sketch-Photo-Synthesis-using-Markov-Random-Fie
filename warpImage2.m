function [warped_img, p1, p2, p3, p4] = warpImage2(shape_mean_scaled, texture_base, triangles, resolution, shape, img, interpolation)
  
  % warpImage2 Summary of this function goes here
  % Detailed explanation goes here

  % Preamble
  n_triangles = size(triangles,1);
  height = size(img, 1);
  width = size(img, 2);

  % Compute gray img
  if (size(img,3) > 1)
      target_gray = rgb2gray(img);   
  else 
      target_gray = img;   
  end

  % Initialize piece-affine affine warp coefficients
  a = zeros(1,6);

  % Initialize warped_image and p1, p2, p3, p4
  warped_img = zeros(resolution(1), resolution(2));
  p1 = zeros(resolution(1), resolution(2));
  p2 = zeros(resolution(1), resolution(2));
  p3 = zeros(resolution(1), resolution(2));
  p4 = zeros(resolution(1), resolution(2));

  for t = 1:n_triangles

    % Compute a = [a1,a2,a3,a4,a5,a6] for each triangle in the shape

    % Coordinates of three vertices of a triangle in reference frame
    U = shape_mean_scaled(triangles(t,:),1);
    V = shape_mean_scaled(triangles(t,:),2);

    % Coordinates of three vertices of a triangle in shape
    X = shape(triangles(t,:),1);
    Y = shape(triangles(t,:),2);

    denominator = (U(2) - U(1)) * (V(3) - V(1)) - (V(2) - V(1)) * (U(3) - U(1));

    a(1) = X(1) + ((V(1) * (U(3) - U(1)) - U(1)*(V(3) - V(1))) * (X(2) - X(1)) + (U(1) * (V(2) - V(1)) - V(1)*(U(2) - U(1))) * (X(3) - X(1))) / denominator;
    a(2) = ((V(3) - V(1)) * (X(2) - X(1)) - (V(2) - V(1)) * (X(3) - X(1))) / denominator;
    a(3) = ((U(2) - U(1)) * (X(3) - X(1)) - (U(3) - U(1)) * (X(2) - X(1))) / denominator;

    a(4) = Y(1) + ((V(1) * (U(3) - U(1)) - U(1) * (V(3) - V(1))) * (Y(2) - Y(1)) + (U(1) * (V(2) - V(1)) - V(1)*(U(2) - U(1))) * (Y(3) - Y(1))) / denominator;
    a(5) = ((V(3) - V(1)) * (Y(2) - Y(1)) - (V(2) - V(1)) * (Y(3) - Y(1))) / denominator;
    a(6) = ((U(2) - U(1)) * (Y(3) - Y(1)) - (U(3) - U(1)) * (Y(2) - Y(1))) / denominator;

    % Warp image into the reference frame through shape base
    max_indexes_img = width * height;

    [v,u] = find(texture_base == t);

    if (~isempty(v) && ~isempty(u))

      indexes_base = v + (u-1) * resolution(1);

      target_pixels_x = a(1) + a(2) .* u + a(3) .* v;
      target_pixels_y = a(4) + a(5) .* u + a(6) .* v;

      indexes_img = round(target_pixels_y) + (round(target_pixels_x)-1) * height;

      % Find pixel coordinates in img
      indexes_img = indexes_img(indexes_img > 0 & indexes_img <= max_indexes_img);

      switch interpolation

        case 'none'

          % Pixel translation
          warped_img(indexes_base) = target_gray(indexes_img);

        case 'bilinear'

          % Bilinear interpolation
          delta_x = target_pixels_x - floor(target_pixels_x);
          delta_y = target_pixels_y - floor(target_pixels_y);
          
          warped_img(indexes_base) =  (1-delta_x) .* (1-delta_y) .* target_gray((floor(target_pixels_y))   + (floor(target_pixels_x)-1) * height) +...
                                       delta_x    .* (1-delta_y) .* target_gray((floor(target_pixels_y))   + (floor(target_pixels_x))   * height) +...
                                      (1-delta_x) .*  delta_y    .* target_gray((floor(target_pixels_y)+1) + (floor(target_pixels_x)-1) * height) +...
                                       delta_x    .*  delta_y    .* target_gray((floor(target_pixels_y)+1) + (floor(target_pixels_x))   * height);

      end
      
      p1(indexes_base) = a(2);
      p2(indexes_base) = a(5);
      p3(indexes_base) = a(3);
      p4(indexes_base) = a(6);

    end
    
  end
  
end

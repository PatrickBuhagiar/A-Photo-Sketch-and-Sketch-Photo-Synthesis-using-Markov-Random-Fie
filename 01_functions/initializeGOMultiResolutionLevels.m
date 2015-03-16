function [dx, dy, p_s] = initializeGOMultiResolutionLevels(param, img, shape_ev, shape_dims, noise_exp)

  % initializeGOMultiResolutionLevels Summary of this function goes here
  % Detailed explanation goes here
  
  % Number of shape eigenvalues
  N_s = length(shape_ev);
  
  % Initialize per level structures
  dx = cell(param.levels, 1);
  dy = cell(param.levels, 1);
  p_s = cell(param.levels, 1);
      
  for i = 1:param.levels

    % Define smoothing parameter sigma
    sigma = i+1; %2^(i-1);

    % If gaussian smoothing, define smoothing kernel and smooth
    if strcmp(param.smoothing, 'gaussian')
      % Kernel
      siz = max(2*sigma+1, 5);
      H  = fspecial('gaussian', [siz siz], sigma);
      % Apply smoothing
      img2 = imfilter(img, H, 'replicate');
    end

    % Compute image gradients
    switch param.gradient

      case 'image'

        if strcmp(param.smoothing, 'longRange')
          % Compute gradients of input image
          [dx{i} dy{i}] = longRangeGradient(img, sigma);
        else
          % Compute gradient
          [dx{i} dy{i}] = gradient(img2);
        end

      case 'rescaled'

        if strcmp(param.smoothing, 'longRange')
          % Compute gradients of input image
          [dx{i} dy{i}] = longRangeGradient(img, sigma);
        else
          % Compute gradient
          [dx{i} dy{i}] = gradient(img2);
        end

    end

    % Compute noise parameter p
    p_s{i} = 1 / (N_s - shape_dims{i}) * sum(shape_ev(shape_dims{i} + 1:end));
    p_s{i} = p_s{i}^noise_exp;
  
  end          
  
end
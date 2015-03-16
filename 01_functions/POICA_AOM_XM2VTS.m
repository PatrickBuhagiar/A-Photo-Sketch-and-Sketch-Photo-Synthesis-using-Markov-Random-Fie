function [final_shape] = POICA_AOM_XM2VTS(param, AAM, image_path)

  % POICA_AOM Summary of this function goes here
  % Detailed explanation goes here
  
  %% Model Parameters
  
  % TEXTURE OPTIONS

  % Levels of multiresolution
  levels = {1, 2, 3, 4};
  param.AAM.levels = cell2mat(levels(4));

  % Interpolation
  interpolation = {'none', 'bilinear'};
  param.AAM.interpolation = cell2mat(interpolation(2));

  % Smoothing
  smoothing = {'gaussian', 'longRange'};
  param.AAM.smoothing = cell2mat(smoothing(1));

  % Gradient Computation
  gradients = {'image', 'rescaled'};
  param.AAM.gradient = cell2mat(gradients(1));
  


  %% Fitting Parameters

  % GENERAL OPTIONS

  % Iterations
  param.n_iterations = 50;

  % Misalignment Offset
  misalignment = {'detector'};
  param.misalignment = cell2mat(misalignment(1));


  % SHAPE OPTIONS

  % Shape Projection
  shape_projection = {true, false};
  param.shape_projection = cell2mat(shape_projection(1));

  % Shape Regularization 
  shape_regularization = {true, false};
  param.shape_regularization = cell2mat(shape_regularization(1));


  % SHAPE & TEXTURE OPTIONS

  % Noise Exponent
  noise_exp = {1, 2/3, 1/2, 1/8};
  param.noise_exp = cell2mat(noise_exp(2));
  
  
  %% Load Image
  
  test_image = imread(image_path);
  

  %% Initialization
 
  % Create video object
  if param.video
    vidObj = VideoWriter(['.' strtok(image_path,'.') '.avi']);
    vidObj.FrameRate = 10;
    vidObj.Quality = 100;
    open(vidObj);
  end
  
  % Define starting multiresolution level
  level = param.AAM.levels;
  change_level_s = true;
  n_iterations_level = param.n_iterations / param.AAM.levels;

  % Test image resolution and number of color channels
  [~, ~, n_channels] = size(test_image); 
  
  % Define input image
  if n_channels == 3
    input_image = double(rgb2gray(test_image))/255;
  else
    input_image = double(test_image)/255;
  end
  
  % Initialize current multiresolution level
  [dx, dy, p_s] = initializeGOMultiResolutionLevels(param.AAM, input_image, AAM.shape_ev, AAM.shape_dims, param.noise_exp);
   

  %% Compute Current Shape

  switch param.misalignment
   
    case 'detector'
        
      % Run Matlab's standard face detector
      faceDetector = vision.CascadeObjectDetector();
      bbox = step(faceDetector, test_image);
      
      % Check if face has been detected
      if isempty(bbox)
        
        face_detected = false;
        
      else
        
        face_detected = true;
        
        % Place shape mean at detected point
        over_scale = 0.9 * bbox(end,3) / AAM.resolution{1}(1);
        current_shape = AAM.shape_mean_level{level} * over_scale - repmat(mean(AAM.shape_mean_level{level} * over_scale, 1), AAM.n_vertices, 1) +  repmat([(bbox(end,1)+bbox(end,3)/2), ((bbox(end,2)+bbox(end,4)/1.65))], AAM.n_vertices, 1) / 2^(level-1);
  
        if param.display
          if n_channels == 3
            boxInserter  = vision.ShapeInserter('BorderColor', 'Custom', 'CustomBorderColor', [0 255 0]);
          else
            boxInserter  = vision.ShapeInserter('BorderColor', 'Custom', 'CustomBorderColor', 0); 
          end
          image_face = step(boxInserter, test_image, uint32(bbox(end,:))); 
          h1 = figure(1);
          imshow(image_face);
          title('Face detection'); 
          drawnow;
          if param.video
            vidObj = saveFrame(vidObj, 5);
          end
          
          if param.save_image
            print(h1, '-dpng', ['.' strtok(image_path,'.') '_fitted0' '.png']);
          end
        end
        
      end
      
  end
 
  if face_detected
    
    if param.verbose
      disp('Face detected');
    end
    
    % Define current_shape2
    current_shape2 = current_shape * 2^(level-1);

    if param.display
      figure(1);
      clf;
      imshow(test_image, []);
      hold on;
      plot(current_shape2(:,1), current_shape2(:,2), 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'MarkerSize', 1.5);
      hold off;    
      drawnow;
      if param.video
        vidObj = saveFrame(vidObj, 5);
      end 
    end
    
    
    %% Fitting
        
    for i = 1:param.n_iterations    
      %% Shape projection
      
      if param.shape_projection && change_level_s
        [current_shape, ~, ~] = computeShapeProjection(current_shape, AAM.shape_mean_level{level}, AAM.shape_pc{level}, AAM.shape_ev, AAM.global_transform{level}, AAM.shape_dims{level}, AAM.triangles, AAM.sorted_triangles{level}, p_s{level}, AAM.n_vertices, param.shape_regularization);
        change_level_s = false;
      end


      %% Warp Input Image and Compute Gradient Orientation
      
      switch param.AAM.gradient

        case 'image'

          % Warp gradients
          hx = warpImage(AAM.shape_mean_scaled{level}, AAM.texture_base{level}, AAM.triangles, AAM.resolution{level}, current_shape2, dx{level}, param.AAM.interpolation);
          hy = warpImage(AAM.shape_mean_scaled{level}, AAM.texture_base{level}, AAM.triangles, AAM.resolution{level}, current_shape2, dy{level}, param.AAM.interpolation);

        case 'rescaled'

          % Warp gradients
          [aux1 p1 p2 p3 p4] = warpImage2(AAM.shape_mean_scaled{level}, AAM.texture_base{level}, AAM.triangles, AAM.resolution{level}, current_shape2, dx{level}, param.AAM.interpolation);
          aux2 = warpImage(AAM.shape_mean_scaled{level}, AAM.texture_base{level}, AAM.triangles, AAM.resolution{level}, current_shape2, dy{level}, param.AAM.interpolation);

          % Rescale
          hx = aux1 .* p1 + aux2 .* p3;
          hy = aux1 .* p2 + aux2 .* p4;

      end

      % Remove magnitude
      ang = angle(hx + 1i * hy);
      ix = cos(ang);
      iy = sin(ang);

      % Remove non-facial pixels
      ix = ix(AAM.texture_mask{level});
      iy = iy(AAM.texture_mask{level});

      im = [ix; iy];
      
      
      %% Gradient Based Correlation Coefficient Cost Funtion

      % - (1) Compute ECC cost function terms -----------------------------

      im = im - AAM.A{level} * (AAM.A{level}' * im);
      u_bold = AAM.G{level}' * im - AAM.G{level}' * AAM.A{level} * (AAM.A{level}' * im);
      u = AAM.t{level}' * im - AAM.t{level}' * AAM.A{level} * (AAM.A{level}' * im);
      den = u - u_bold' / AAM.Q{level} * AAM.v_bold{level};


      % - (2) Compute parameter lambda ------------------------------------

      if u >  u_bold' / AAM.Q{level} * AAM.v_bold{level}
        lambda_shape = AAM.num{level} / den;
      else
        if u <=  u_bold' / AAM.Q{level} * AAM.v_bold{level}
          lambda_shape1 = sqrt((AAM.v_bold{level}' / AAM.Q{level} * AAM.v_bold{level}) / (u_bold' / AAM.Q{level} * u_bold));
          lambda_shape2 = (u_bold' / AAM.Q{level} * AAM.v_bold{level} - u) / (u_bold' / AAM.Q{level} * u_bold);
          lambda_shape = max(lambda_shape1, lambda_shape2);
        else
          disp('Error in computation of lambda');
          break;
        end
      end


      % - (3) Compute parameter updates -----------------------------------

      % Compute delta SDs
      SD_delta = lambda_shape * u_bold - AAM.v_bold{level};

      % Compute parameter updates
      delta_rhos_alphas = AAM.Q{level} \ SD_delta;


      %% Shape Compositional Update

      % Shape composition
      [current_shape, ~, ~] = computeWarpUpdate(current_shape, delta_rhos_alphas, AAM.shape_mean_level{level}, AAM.shape_pc{level}, AAM.shape_ev, AAM.global_transform{level}, AAM.n_transformations, AAM.shape_dims{level}, AAM.triangles, AAM.sorted_triangles{level}, p_s{level}, AAM.n_vertices, param.shape_regularization);
      
      % Define current_shape2
      current_shape2 = current_shape * 2^(level-1);


      %% Display

      if param.display
        % Show fitting evolution
        figure(1);
        clf;
        imshow(test_image, []);
        hold on;
        trimesh(AAM.triangles, current_shape2(:,1), current_shape2(:,2), 'g', 'MarkerEdgeColor', 'g', 'MarkerSize', 2);
        plot(current_shape2(:,1), current_shape2(:,2), 'o', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'MarkerSize', 2);
        hold off;
        drawnow;
        if param.video
          vidObj = saveFrame(vidObj, 1);
        end
      end
      

      %% Determine Level Convergence

      if i >  (param.AAM.levels+1 - level) * n_iterations_level
        % Upgrade level
        level = level - 1;
        change_level_s = true;
        current_shape = current_shape * 2;
      end


    end


    %% Loop Postscript

    if param.display
      % Show ground truth and last iteration shapes
      h2 = figure(2);
      clf;
      imshow(test_image, []);
      hold on;
      trimesh(AAM.triangles, current_shape2(:,1), current_shape2(:,2), 'g', 'MarkerEdgeColor', 'g', 'MarkerSize', 2);
      plot(current_shape2(:,1), current_shape2(:,2), 'o', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'MarkerSize', 2);
      title('Final fitting');
      hold off;
      drawnow;
      if param.video
        vidObj = saveFrame(vidObj, 5);
      end
    end

    final_shape = current_shape2;
    
    
  else
    
    if param.verbose
      disp('No face detected');
    end
  
    final_shape = zeros(AAM.n_vertices, 2);  
    
  end
  
  
  %% Postscript

  if param.video
    close(vidObj);  
  end
  
  if param.save_image
    print(h2, '-dpng', ['.' strtok(image_path,'.') '_fittedX' '.png']);
  end


end
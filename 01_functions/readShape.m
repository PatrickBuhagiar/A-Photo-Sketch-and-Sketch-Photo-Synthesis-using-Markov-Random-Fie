function [shape] = readShape(shape_path, n_vertices)

  % readShape Summary of this function goes here
  % Detailed explanation goes here
  
  [~, remain] = strtok(shape_path,'.');

  switch remain

    case '.asf'

      % Load shape
      shape = dlmread(shape_path, '\t', [16 2 16+n_vertices-1 3]);
      shape(:,1) = shape(:,1);
      shape(:,2) = shape(:,2);

    case '.pts'

      % Skip unrellevant fields
      fid = fopen(shape_path);
      tline = fgetl(fid);
      start = 1; 
      while ~strcmp(tline, '{')     
          start = start + 1;
          tline = fgetl(fid);      
      end
      fclose(fid);

      % Load shape
      shape =  dlmread(shape_path, ' ', [start 0 start+n_vertices-1 1]);
      
    case '.pts3'

      % Skip unrellevant fields
      fid = fopen(shape_path);
      tline = fgetl(fid);
      start = 1; 
      while ~strcmp(tline, '{')     
          start = start + 1;
          tline = fgetl(fid);      
      end
      fclose(fid);

      % Load shape
      shape =  dlmread(shape_path, ' ', [start 0 start+n_vertices-1 2]);

  end
    
end
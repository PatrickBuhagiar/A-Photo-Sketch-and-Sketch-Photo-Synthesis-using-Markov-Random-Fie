function [vidObj] = saveFrame(vidObj, n_copies)

  % saveFrame Summary of this function goes here
  % Detailed explanation goes here
  
  for j = 1:n_copies
    frame = getframe;
    writeVideo(vidObj, frame);
  end

end
function [d_current, d_neighbour]= get_overlap(current_patch_index, current_patch, neighbour_patch_index, neighbour_patch, overlap, full_image_size)
    [y,x] = size(current_patch);
    z = full_image_size(2);
    Xend = round(x:(x - overlap):z);
    width = length(Xend);
    %if up
    if (current_patch_index == neighbour_patch_index + width)
        d_current = current_patch(1:overlap, :)';
        d_neighbour = neighbour_patch(y-overlap+1:1:y, :)';
    %if down
    elseif (current_patch_index == neighbour_patch_index - width)
        d_current = current_patch(y-overlap+1:1:y, :)';
        d_neighbour = neighbour_patch(1:overlap, :)';
    %if left
    elseif (current_patch_index == neighbour_patch_index + 1)
        d_current = current_patch(:, 1:1:overlap);
        d_neighbour = neighbour_patch(:,(x-overlap+1):1:x);
    %if right
    elseif (current_patch_index == neighbour_patch_index - 1)
        d_current = current_patch(:,(x-overlap+1):1:x);
        d_neighbour = neighbour_patch(:, 1:1:overlap);
    end


end
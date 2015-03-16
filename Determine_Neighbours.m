function neighbours = Determine_Neighbours(n_patches, patch_size, overlap_size, full_image_size)
    %n_patches: total number of patches
    %full_image_size: [x,y]
    x = full_image_size(2);
    Xend = round(patch_size:(patch_size - overlap_size):x);
    width = length(Xend);
    
    corner = [1, width, (n_patches - width -1), n_patches];
    left_edge = [(width+1):width:(n_patches -(width*2)+1)];
    right_edge = [(width*2):width:(n_patches-width)];
    
    neighbours = cell(1, n_patches);
    
    for i=1:n_patches,     
        %% corners
        %if top left corner
        if (i==1)
            temp = [i+1 i+width];
        %if top right corner
        elseif (i==width)
            temp = [i-1 i+width];
        %if bottom left corner
        elseif (i==(n_patches - width -1))
            temp = [i-width i+1];  
        %if bottom right corner
        elseif (i==n_patches)
            temp = [i-width i-1];
        %% top and bootom edges
        %if top
        elseif (i<width)
            temp = [i-1 i+1 i+width];
        %if bottom
        elseif (i>(n_patches-width))&&(i<n_patches)
            temp = [i-width i-1 i+1];
        %% Edges
        %if left edge
        elseif (is(left_edge, i))
            temp = [i-width i+1 i+width];
        %if right edge
        elseif (is(right_edge, i))
            temp = [i-width i-1 i+width];
        else 
            temp = [i-width i-1 i+1 i+width];
        end
        neighbours{i}=temp;
    end    
end

function b = is(array, index)
    if any(index==array)
        b = true;
    else 
        b = false;
    end
end
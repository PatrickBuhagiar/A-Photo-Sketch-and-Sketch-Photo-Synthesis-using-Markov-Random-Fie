function imOut = GeneratePsuedoPhoto(Patches, ImageSize, overlapSize)
    
    %Initialise Output
    x = ImageSize(1);
    y = ImageSize(2);
    imOut = zeros(ImageSize);
    patchSize = size(Patches, 1);
    
    %Resize image to fit patches
    Xstart = 1:(patchSize - overlapSize):x;
    Ystart = 1:(patchSize - overlapSize):y;   
    Xend = round(patchSize:(patchSize - overlapSize):x);
    Yend = round(patchSize:(patchSize - overlapSize):y); 
    imOut = SubSample(imOut,Xend(end),Yend(end));
    
    %start and end bounds for overlaps
    XoverlapStart = (patchSize - overlapSize + 1):(patchSize - overlapSize):x;
    XoverlapEnd = XoverlapStart + 4;
    YoverlapStart = (patchSize - overlapSize + 1):(patchSize - overlapSize):x;
    YoverlapEnd = XoverlapStart + 4;
    
    
    PatchesCount = [size(Xend, 2) size(Yend,2)]; 
    c = 1;
    for a=1:PatchesCount(1),
        for b=1:PatchesCount(2),
            if b==1 && a==1 %first patch. Not stitching required
                imOut(Xstart(a):1:Xend(a),Ystart(b):1:Yend(b)) = Patches(:,:,c);
            elseif b~=1 && a==1 %first row of images. stitching required on one side.
                %Temporary holder for patch 
                canvas = Patches(:,:,c);  
                %Calculate Overlap (left side)
                overlap = (imOut(Xstart(a):1:Xend(a),Ystart(b):1:Ystart(b)+(overlapSize-1)) + canvas(:,1:1:(overlapSize)))./2;
                %place patch
                imOut(Xstart(a):1:Xend(a),Ystart(b):1:Yend(b)) = Patches(:,:,c);
                %replace overlap area with calculated overlap
                imOut(Xstart(a):1:Xend(a), Ystart(b):1:Ystart(b)+(overlapSize-1)) = overlap; 
            else
                %Temporary holder for patch 
                canvas = Patches(:,:,c);  
                %Calculate Overlap (left side)
                overlap = (imOut(Xstart(a):1:Xend(a),Ystart(b):1:Ystart(b)+(overlapSize-1)) + canvas(:,1:1:(overlapSize)))./2;
                %Calculate Overlap (top side)
                overlap2 = (imOut(Xstart(a):1:Xstart(a)+(overlapSize-1),Ystart(b):1:Yend(b)) + canvas(1:1:(overlapSize),:))./2;
                %place patch
                imOut(Xstart(a):1:Xend(a),Ystart(b):1:Yend(b)) = Patches(:,:,c);
                %replace overlap area with calculated overlap
                imOut(Xstart(a):1:Xend(a), Ystart(b):1:Ystart(b)+(overlapSize-1)) = overlap; 
                imOut(Xstart(a):1:Xstart(a)+(overlapSize-1),Ystart(b):1:Yend(b)) = overlap2; 
            end
            %imOut(Xstart(a):1:Xend(a),Ystart(b):1:Yend(b)) = (imOut(Xstart(a):1:Xend(a),Ystart(b):1:Yend(b)) + Patches(:,:,c));  
  
        c = c+1;
        end
    end
    
%     for a=1:(PatchesCount(1)-1),
%         for b=1:(PatchesCount(2)-1),
%             imOut(XoverlapStart(a):1:XoverlapEnd(a),YoverlapStart(b):1:YoverlapEnd(b)) = imOut(XoverlapStart(a):1:XoverlapEnd(a),YoverlapStart(b):1:YoverlapEnd(b))./2;
%         end
%     end
%     
end
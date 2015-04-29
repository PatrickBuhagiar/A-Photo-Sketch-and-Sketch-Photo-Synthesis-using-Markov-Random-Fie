% function face_rec
%
% Train function performing facial identification
%
% Inputs:
%   train_data - data to be used for training (e.g. if doing face-to-sketch transformation, need array containing sketches on which to train FR)
%   fname - Name of face recognition method (e.g. 'Eigenfaces')
%   mode - mode to use for facial identifier ('training'/'testing')
%
% Outputs:
%   outData - results, containing best matches in descending order (i.e. best match is first, next best match is second etc.)

function outData = face_rec(train_data, fname, mode)

if strcmpi(fname, 'Eigenfaces')
    outData = Eigenfaces5(train_data, mode); % Model of transformation using eigentransformation
else
    error('This Face Recognition Method is not supported yet');
end

% function FR_train_data = face_rec(train_data_FR, FR_fn, mode)
% 
% if strcmpi(FR_fn, 'Eigenfaces')
%     FR_train_data = Eigenfaces5(train_data_FR, mode); % Model of transformation using eigentransformation
% else
%     error('This Face Recognition Method is not supported yet');
% end

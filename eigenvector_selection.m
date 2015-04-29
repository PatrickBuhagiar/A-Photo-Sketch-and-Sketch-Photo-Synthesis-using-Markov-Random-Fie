function [V,D] = eigenvector_selection(Evector,Evalues, T)
%tic
% Sort the Eigenvalues in descending order
[Evalues, I] = sort(Evalues,'descend');

% Initialize the accumulative eigenvalues
s = 0;
for n = 1:size(Evalues,1)
    % Derive the accumulative variance
    s = s + Evalues(n);
    
    % Calculate the overall energy within the current n
    ratio = s/sum(Evalues,1);
    
    if ratio > T
        % Truncate the other indices
        I = I(1:n);
        % Return the vale
        break;
    end
end

% Derive the eigenvectors selected
V = Evector(:,I);

% Derive the selected eigenvalues
D = Evalues(size(Evalues,1)-I+1);%(1:n);
%toc
%% New eigenvector selection
% 
% V=[];
% D=[];
% for i=1:size(Evector,2)
%     if(Evalues(i)>1e-4)
%         V=[V Evector(:,i)];
%         D=[D Evalues(i)];
%     end
% end

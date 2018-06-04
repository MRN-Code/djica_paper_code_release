function [W,runtimes] = dpca(localVectors,comp_number,v_)
X = localVectors{1}; % assuming voxel x time matrix here
colsX = size(X,2);
fprintf('\t\t\t\tStacking and reducing...\n');
fprintf('\t\t\t\tnodes done: ');
for site = 2:length(localVectors)
    tic
    X = [X localVectors{site}]; % Stack vectors in component dimension
    [H,~] = eig(X'*X);
    X = X*H(:,end:-1:(end-colsX+1));    % Reduce back to 'colsX' columns
    if mod(site,1) == 0
        fprintf('/%i',site);
    end
end
fprintf('\n');

W = X(:,1:comp_number); % Final number of global PCs to keep
fprintf('\t\t\t\tNormalizing...\n');
n = zeros(1,size(W,2));
for kk = 1:size(W,2)
    n(kk) = norm(W(:,kk));  % Norm of each column
end
W = W*diag(n.^(-1));

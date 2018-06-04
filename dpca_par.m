function W = myGlobalPCA_P(localVectors,comp_number)
numSites=length(localVectors);
clusterSize=2;
numCluster=floor(numSites/clusterSize);
fprintf('Number of Sites Remaining %i, Size of Clusters %i, Number of Clusters %i\n', numSites, clusterSize, numCluster);
clusterW = cell(1,numCluster);
if numCluster > 1
    jump = brad_compute_jump(numCluster, 256, 2);
    for ii = 1:jump:numCluster
        parfor cluster = ii:min(ii+(jump-1), numCluster);
            serialVect = localVectors;
            beginningSite = (cluster-1)*clusterSize+1;
            endSite = (cluster)*clusterSize;
            if cluster==numCluster
                endSite = numSites;
            end
            lV=serialVect(beginningSite:endSite);
            clusterW{cluster} = myGlobalPCA_P(lV, comp_number);
        end
    end
    W = myGlobalPCA_P(clusterW, comp_number);
else
X = localVectors{1}; % assuming voxel x time matrix here
colsX = size(X,2);

fprintf('\t\t\t\tStacking and reducing...\n');
fprintf('\t\t\t\tnodes done: ');
for site = 2:length(localVectors)
    X = [X localVectors{site}]; % Stack vectors in component dimension
    %spmd
%	C = codistributed(X'*X);
	[H,~] = eig(X'*X);
    %end
    %H = H{:};
    X = X*H(:,end:-1:(end-colsX+1));    % Reduce back to 'colsX' columns
    if mod(site,1) == 0
        fprintf('/%i',site);
    end
end
fprintf('\n');
W=X;

end
end

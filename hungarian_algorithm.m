function varargout = hungarian_algorithm(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reassign_sources - reassigns sources using the hungarian algorithm
%
%   usage >> [src1, src2,...,srcN] =
%   hungarian_algorithm(src1,Corr1,...,srcN,CorrN)
%
%   where srcN is the matrix of sources nand CorrN is the correlation
%   matrix with the ground truth
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% easy arg access %%%%%%%%%%%%%%%%%%%%
infoN = @(index) eval(['src' int2str(index)]);
setInfoN = @(name,index,value) eval([name int2str(index) '=' value]);
%

if mod(length(varargin)/2,1) ~= 0
    fprintf('Please enter a correlation matrix for each dataset\n');
end

datasets = cell(1,length(varargin)/2);
corrs = cell(1,length(varargin)/2);

datasets = varargin(1:length(datasets));
corrs = varargin(length(datasets)+1:end);
j = 1;
for i = 1:length(datasets)
    ASSIGN = munkres(-corrs{j});
    tmp = datasets{i}(ASSIGN,:,:);
    varargout{i} = tmp;
    if length(corrs) > j 
        j = j+1;
    end
end

end

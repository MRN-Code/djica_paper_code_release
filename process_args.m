function [keys,vals,datasets,knl,vnl] = process_args(argincell,dataflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Processes arguments given pairwise, i.e.
%	{keyword,value,keyword,value,...,keyword,value}
%
%also handles pairwise arguments including datasets
%	
%usage
%	>> [keys,vals] = process_args({key1,value1,key2,value2,...,keyN,valueN})
%	returns keys and values, treating datasets as values corresponding to keys
%	N must be even!
%	>> [keys,vals,datasets] = process_args({key1,value1,key2,value2,...,keyN,valueN},TRUE)
%	returns keys, values, and datasets separately
%	N can be odd!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('dataflag','var')
    dataflag = 1;
end
    %% Arg Processing
if dataflag
    % Count and collect the passed datasets
    data_count = 0;
    indices = zeros(1,length(argincell));
    for i = 1:length(argincell)
        if (isnumeric(argincell{i}) && size(argincell{i},1) > 1 && size(argincell{i},2) > 1)
            data_count = data_count + 1;
            indices(i) = i;
        end
    end
    indx = indices(indices ~= 0);

    % Place each dataset into a cell
    datasets = cell(1,data_count);
    for i = 1:data_count
        datasets{i} = argincell{indx(i)};
    end
    %end data processing

    % Place keywords and values into cells
    tmp = 1:length(argincell);
    indices = tmp(tmp ~= indices);
    key_count = floor((length(argincell) - data_count)/2);
    if mod(key_count,1) ~= 0
        fprintf('distpca():-keywords-please put keywords and values in pairs');
        return;
    end
    keys = cell(1,key_count);
    vals = cell(1,key_count);
    j = 1;
    for i = 1:2:key_count*2
        if i == 0
            break;
        end
        keys{j} = argincell{indices(i)};
        vals{j} = argincell{indices(i+1)};
        j = j + 1;
    end
else
    keys = cell(1,length(argincell)/2);
    vals = cell(1,length(argincell)/2);
    j = 1;
    for i = 1:2:length(argincell)
        keys{j} = argincell{i};
        vals{j} = argincell{i+1};
        j = j +1;
    end
end
i = 1; j = 1;
knl = keys;
vnl = vals;
while (i <= length(keys) && j <= length(vals))
    if ischar(keys{i})
        knl{i} = keys{i};
        keys{i} = lower(keys{i});
    end
    if ischar(vals{j})
        vnl{i} = vals{i};
        vals{j} = lower(vals{j});
    end
    i = i + 1;
   j = j + 1;
end
%end key and value processing 
end

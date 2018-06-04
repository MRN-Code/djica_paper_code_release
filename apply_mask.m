function [varargout] = apply_mask(varargin)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% applies a given mask to a dataset, or calculates and applies mask to a given dataset
	%
	%	usage
	%		>> datasets = apply_mask(datasets,'mask',mask);
	%		where datasets is a cell containing simtb datasets
	%		and mask is a given mask for those datasets
	%		this applies the mask and returns a cell containing the masked datasets
	%		or
	%		>> [mask, datasets] = brad_apply_mask(datasets);
	%		this returns the mask and a cell containing the masked datasets
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% args
    [keys, vals,datasets] = process_args(varargin,1);
    flag_demean = 0;
    for i = 1: length(keys)
        Key = lower(keys{i});
        Val= vals{i};
        if ischar(Val)
            Val = lower(Val);
        end
        
        if strcmp(Key,'mask')
            if ~isnumeric(Val) || (size(Val,1) > 1 && size(Val,2) > 1)
                error('dICA:apply_mask:ValueInputError','Mask must be numerical vector');
            else
                mask = Val;
                varargout = cell(1,length(datasets));
                len = 0;
            end
        elseif strcmp(Key,'data')
            if ~isnumeric(Val) || (size(Val,1) < 1 && size(Val,2) < 1)
                error('dICA:apply_mask:ValueInputError','data must be numerical matrix');
            else
                if ~iscell(Val)
                    Val = {Val};
                end
                datasets = Val;
                varargout = cell(1,length(datasets));
                len = 0;
            end
        elseif strcmp(Key,'demean')
            flag_demean = Val;
        end
        
    end

	%%%%%%%%%%%%%%%%%% MASKING %%%%%%%%%%%%%%%%%%%%%%

	%calculate individual masks for subjects
    if ~exist('mask','var')
    	mask = brad_calc_mask(datasets{:});
        varargout{1} = mask;
        len = 1;
    end
 
	%apply mask
    for i = 1:length(datasets)
        try
            if flag_demean
                varargout{len + i} = detrend(datasets{i}(mask==1,:),0);
            else
                varargout{len + i} = datasets{i}(mask==1,:);
            end
        catch err
            if ~flag_demean
                varargout{len + i} = datasets{i}(:,mask==1);
            else
                varargout{len + i} = detrend(datasets{i}(:,mask==1),0);
            end
        end
    end
    %
    
end

function [mask,subjmasks,counts] = compute_mask...
    (varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc_mask - calculates the voxel-mask between multiple datasets,        
% multiplying together individual data-masks element-wise to create the
% full datamask. 
%
%   usage  >>mask = compute_mask(data1,data2,data3,...,data15);
%
%   use flag 'mask' to pass own mask 
%
%   input data - data1,...,data15 - takes datasets of arbitrary size and
%   then calculates their masks individually, multiplying element-wise to
%   obtain a mask between all datasets. 
%
%   output data - mask - mask is the unapplied mask for all datasets
%   entered. Entering one dataset will calculate the mask for only that
%   dataset.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Edit history %%%%%%%%%%

%   created on Friday June 27th
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% command-line variables
    flag_transpose = 1;
    maskType = 'mean';
    flag_dim = 'row'; % reduce in the column dimension
    
    %% argument processing 
    [keys,vals,datasets] = process_args(varargin,1);
    
    argK = {'maskType',     @(v)ischar(v);...
            'flag_transpose', @(v) (v == 1 || v == 0);...
            'flag_dim', @(v) ischar(v) && (strcmp(v,'col') || strcmp(v,'row'));...
           };
    
    %%% Create Map container for access by keywords
    V=cell(1,size(argK,1));
    K=cell(1,size(argK,1));
    for i = 1:size(argK,1)
        V{i} = i;
        K{i} = lower(argK{i,1}); %needs to be lowered for comparison
    end
    argMap = containers.Map(K,V);
    %%% Iterate through keywords
    for i = 1:length(keys)
        if (isKey(argMap,keys{i}) ... % is the key in our map
                && argK{argMap(keys{i}),2}(vals{i}) ... %  is it logical?
                )
             eval([argK{argMap(keys{i}),1} ' = vals{i};']); %set variable
        end
    end 
    %% Post - Processing  
    [m n] = cellfun(@size,datasets);
    n = min(n);
    m = min(m);
    rdim = 1;
    if strcmp(flag_dim,'col')
        rdim = 2;
    end
    %% initialize mask as ones
    mask = ones(size(datasets{1},rdim),1);
    counts = zeros(1,length(datasets));
    subjmasks = cell(size(datasets));

    %% Mean Mask - NaN Mask - Circle Mask (default)
    if strcmp(maskType,'mean')
        %%%% multiple mean mask calculation loop %%%%
        for i = 1:length(datasets) %all arguments are assumed to be datasets
            mask_n = ones(size(datasets{i}, 2), 1); %initialize for nth dataset

            %%% single mask calculation loop
            for nM = 1%:size(datasets{i}, 1)
                tmp = datasets{i}(nM, :)';
                mask_n = mask_n.*(tmp >= mean(tmp(:)));
            end% end single mask calculation loop
            
            %update the mask for multiple datasets
            counts(i) = sum(~mask_n);
            subjmasks{i} = mask_n;
            mask = mask.*mask_n;

        end%end multiple mask calculation loop
        
%         mask = mask.*mask_n;

        
    elseif strcmp(maskType,'nan')
        %%%% multiple nan mask calculation loop %%%%

        for i = 1:length(datasets)
            mask_n = ones(size(datasets{i}, 2), 1); %initialize for nth dataset

            %%% single mask calculation loop
            for nM = 1:size(datasets{i}, 1)
                tmp = datasets{i}(nM, :)';
                mask_n = mask_n.*(~isnan(tmp));
               
            end% end single mask calculation loop
            counts(i) = sum(~mask_n);
            %update the mask for multiple datasets
            mask = mask.*mask_n;
        end
    elseif strcmp(maskType,'circ')
        d = n;
        if strcmp(flag_dim, 'row')
            d = m;
        end
        nV = sqrt(max(d));
        arg1 = linspace(-1,1,nV);
        [x,y] = meshgrid(arg1,arg1);
        r = sqrt(x.^2 + y.^2);

        mask = ones(1,nV*nV);
        mask(r>1) = 0;
        if flag_transpose
            mask = mask';
        end
    end
    
end

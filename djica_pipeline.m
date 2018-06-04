function varargout = djica_pipeline(datasets, varargin)
    %% djica_pipeline
    %   This function combines all of the steps required in the distribute
    %   joint ICA pipepline for simulated and real data. The scrThe script uses the MATLAB
    %   parallel computing toolbox if available, but defaults to running in serial
    %
    %   USAGE: >>djica_pipeline(datasets, 'key1',val1,'key2',val2,...,'keyN',valN);
    %       where datasets is a cell containing individual subject data 
    % 		'keyN' is the Nth keyword, and valN is the Nth value
    %       associate with that proceeding keyword.
    %
    %    | KEYWORD ; ACCEPTED VALUES ; DEFAULT VALUE ; DESCRIPTION
    %    | %%% Data Information
    %    | subjMat ; Numeric Vector ; [4] ; Number of Subjects
    %    | siteMat ; Numeric Vector ; [2] ; Number of Sites
    %    | 
    %    | %%% ICA and PCA components
    %    | nc ; numeric ; 50 ; Number of Components
    %    | k ; numeric ; 100 ; Dimension for Local PCA
    %    | 
    %    |%%% flags
    %    | flag_local_pca ; 0 or 1; 0; run local pca 
    %    | flag_dpca 	; 0 or 1; 0; run dpca
    %    | flag_djica 	; 0 or 1; 0; run djica
    %    | flag_isi	; 0 or 1; 0; compute the isi - only if djica is run
    %    | flag_verbose ; 0 or 1; 0; verbose display flag
    %    | flag_mask ; 0 or 1; 1; perform spatial map masking
    %    | flag_demean ; 0 or 1; 1; perform distributed demeaning
    %    | flag_timer ; 0 or 1; 1; record times of certain parts of analysis
    %    | flag_normalize ; 0 or 1; 1; perform normalization
    %    |    
    %
    %   The full running of the function requires several external scripts
    %   to function correctly, which are included in the code release.
    %
    %   This program is dependent on MRN's ICA toolbox, and the following
    %   scripts:
    %       process_args.m
    %	    compute_mask.m
    %       apply_mask.m
    %       hungarian_algorithm.m
    %       fill_mask.m
    % 	    compute_isi.m

    
    %% Input Preprocessing
 
    %%%% DEFAULT values for keyword usage on command line
    %%% Experiment Info Defaults
    DEFAULT_subjMat = 16;              % Vector of numbers - must be the same size as siteMat. The program runs through pairs of subject versus site coordinates and records a run for each pair.
    DEFAULT_siteMat = 2;                % Vector of numbers - must be the same size as subjMat.
    DEFAULT_numRuns = 1;                % Number of times to repeat runs, useful for robustness testing using different simulated data, etc.
    
    %%% ICA and PCA component Defaults
    DEFAULT_nc = 20;                    % 20 components for simulated case
    DEFAULT_k = 100;                    %  k is the number of components for the first PCA reductio
    
    DEFAULT_meanDim = 2;		% dimension for demeaning, by default, demean from columns
    
    %%% DEFAULT flags - used for flagging certain areas of code for usage
    DEFAULT_flag_local_pca = 1;				% yes compute local pca
    DEFAULT_flag_dpca = 1;				% yes compute dpca
    DEFAULT_flag_djica = 1;				% yes compute djica
    DEFAULT_flag_isi = 1;				% yes compute isi
    DEFAULT_flag_mask = 1;                              % yes mask voxels outside head
    DEFAULT_flag_demean = 1;                            % yes demean data
    DEFAULT_flag_verbose = 0;                           % yes verbose
    DEFAULT_flag_normalize = 1;                         % yes normalize spatial maps
    DEFAULT_flag_rand_dist = 0;				% do not randomly distribute subjects by default	
    DEFAULT_flag_zscore = 1;				% do compute zscores of maps
    DEFAULT_load_mask = 0;

    %%% DEFAULT save directory
    % this program uses the following file structure for saving experiment
    % results (variable names put in angle brackets <> ):
    % <projDir>/results/<dateDir>/<timeDir>
    %projDir, testPath, dateDir, dataName, MODEL_NAME, timeDir
    DEFAULT_projDir = '.';     % project directory
     
    % Other variables 
    DEFAULT_varyOn       = 'STOP';	% the name of the variable to varyon, this variable can then be entered as a vector,
					% and repeated runs will be computed for the different values
    DEFAULT_distPdf = 'normal';				% pdf used when randomly distributing data
    DEFAULT_nfRuns = 0;					% number of runs to skip in a multi-run setting
    DEFAULT_seed         = 14159;	% default seed
    
    %%%%%%%%%%%%% Pdf Parameters
    DEFAULT_rMean 	 = 128;
    DEFAULT_rStd 	 = 128;

    %%%% Variable Values -- Variables are first initialized to defaults
    %%%% before keyword processing
    %%% Experiment Variables
    subjMat = DEFAULT_subjMat;
    siteMat = DEFAULT_siteMat;
    numRuns = DEFAULT_numRuns;
    distPdf = DEFAULT_distPdf;
    nfRuns = DEFAULT_nfRuns;
    load_mask = DEFAULT_load_mask; 
    %%% ICA and PCA components
    nc = DEFAULT_nc;
    k = DEFAULT_k;
    
    meanDim = DEFAULT_meanDim;
    
    %%% flags
    flag_local_pca = DEFAULT_flag_local_pca;
    flag_dpca = DEFAULT_flag_dpca;
    flag_djica = DEFAULT_flag_djica;
    flag_isi = DEFAULT_flag_isi;
    flag_verbose = DEFAULT_flag_verbose;
    flag_mask = DEFAULT_flag_mask;
    flag_demean = DEFAULT_flag_demean;
    flag_normalize = DEFAULT_flag_normalize;
    flag_rand_dist = DEFAULT_flag_rand_dist;
    flag_zscore = DEFAULT_flag_zscore;
    %%% Save Directory
    projDir = DEFAULT_projDir;
     
    varyOn = DEFAULT_varyOn;
    seed = DEFAULT_seed;
    
    rMean = DEFAULT_rMean;
    rStd = DEFAULT_rStd;

    %% Input Processing 
    [keys vals] = process_args(varargin,0);        % processes arguments given in pairs of keywords and values, e.g. djica_RealRun('key1',val1,'key2',val2,...,'keyN',valN). 
                                                        % type 'help process_args' or see ./process_args.m for more information
    % clear varargin for workspace cleanliness
    clear varargin;
    
    %%% Logicals and IDs for keys - this is a paradigm for assigned and checking the
    %%% logic on user-input variables which keeps us from having huge
    %%% nested if-then blocks and repetitive code. Basically, argK is a cell where each
    %%% row contains information for one variable, the first column is the
    %%% variable name as given in the code, and the second column is the
    %%% logical function which checks the variable's value.
    argK = {'subjMat',         @(v)isnumeric(v);... 
             'siteMat',         @(v)isnumeric(v);...
             'numRuns',         @(v)isnumeric(v);...
             'flag_verbose',    @(v) (v == 1 || v == 0);...
             'flag_mask',       @(v) (v == 1 || v == 0);...
             'flag_demean',     @(v) (v == 1 || v == 0);...
             'flag_normalize',  @(v) (v == 1 || v == 0);...
             'flag_local_pca',  @(v) (v == 1 || v == 0);...
             'flag_dpca',  @(v) (v == 1 || v == 0);...
             'flag_djica',  @(v) (v == 1 || v == 0);...
             'flag_isi',  @(v) (v == 1 || v == 0);...
             'flag_zscore',  @(v) (v == 1 || v == 0);...
             'nc',              @(v)isnumeric(v);...
             'k',              @(v)isnumeric(v);...
	     'nfRuns',  	@(v)isnumeric(v);...
            'distPdf', @(v)ischar(v);...
             'meanDim',         @(v) isnumeric(v) & (v == 1 || v == 2);...
             'varyOn',          @(v) ischar(v);...
             'load_mask',          @(v) ischar(v);...
             'seed',            @(v) isnumeric(v);...
	     'rMean', 		@(v) isnumeric(v);...
	     'rStd',		@(v) isnumeric(v);...
             };
         
    %%% Create Map container for access by the keywords the user inputs.
    %%% The Map data structure allows us to do this access via the keywords
    %%% and keeps us from using unnecessary strcmp commands and huge nestd
    %%% if-then blocks. Using the map structure also allows us to easily
    %%% assign the values associated with keywords in the map to variables
    %%% in our workspace.
    V=cell(1,size(argK,1));
    K=cell(1,size(argK,1));
    for i = 1:size(argK,1)
        V{i} = i;
        K{i} = lower(argK{i,1}); % user input keywords get lowered by process_args, so we need to lower keys for the map
    end
    argMap = containers.Map(K,V);

    %%% Iterate through user-input keywords
    for i = 1:length(keys)
        if strcmp(keys{i},'dataname') % a one-time check to make sure the data name is all-caps. Just a preference on my end.
            vals{i} = upper(vals{i});
        end
        if (isKey(argMap,keys{i}) ... % if the user-input keyword matches any of our accepted keys 
                && argK{argMap(keys{i}),2}(vals{i}) ... %  and it passes the logic of that key
                )
             eval([argK{argMap(keys{i}),1} ' = vals{i};']); % then set the corresponding variable to the user-input value corresponding with that key
        end
    end
    if ~strcmp(varyOn,'STOP')
        varyS = eval(argK{argMap(varyOn),1});
        varyOn = argK{argMap(varyOn),1};
    else
        varyS = 1;
    end

    %%% Clear all of the processing variables for workspace cleanliness
    clearvars DEFAULT* V K argMap argK keys vals 
    %% Input Post-Processing
    %%% Post processing steps include the following steps...
    % 1. initialize experiment saving directory
    % 2. initialize site distribution matrix
    % 3. set correct number of components based on hybrid-data usage
    %%% Parent Directory 
    if ispc
         parent_directory = sprintf('%s\results\',projDir);
    else
        parent_directory = sprintf('%s/results/',projDir);
    end
    [VOX TIME] = cellfun(@size, datasets, 'UniformOutput', 1); 

    if ~exist(parent_directory,'dir')
            mkdir(parent_directory);
    end
    
    for v = 1:numel(varyS)
        randS = RandStream('mt19937ar','Seed',seed);
        rng(seed,'twister');
	for r = 1:nfRuns
		dummy = randperm(randS,length(datasets));
	end
    
        for r = (1+nfRuns):numRuns
            
	    % Shuffle data each run
	    datasets = datasets(randperm(randS, length(datasets)));
            
            if ~strcmp(varyOn,'STOP')
                if iscell(varyS)
                    eval(sprintf([varyOn '=[' num2str(varyS{v}) ']']));  
                else
                    eval(sprintf([varyOn '=[' num2str(varyS(v)) ']']));   
                end
            end

            siteDist = subjMat./siteMat;
            if flag_rand_dist
		% Randomly distribute subjects across sites
                siteDist = cell(1,length(subjMat));
		for s = 1:length(siteDist)
			if strcmp(distPdf,'normal')
				vals = ceil(random('Normal',rMean,rStd,[1 100]));
			elseif strcmp(distPdf,'uniform')
				vals = ceil(random('Uniform',rMean,rStd,[1 100]));
			elseif strcmp(distPdf,'exponential')
				vals = ceil(random('Exponential',rMean,[1,100]));
			else
				vals = ceil(random('Normal',rMean,rStd,[1 100]));
			end
			vals = vals(vals >= 4);
			cs = cumsum(vals);
			last = find(diff(cs > subjMat(s)));
			vals = vals(1:last);
			add = floor((subjMat(s)-sum(vals))/length(vals));
			vals = vals + add*ones(size(vals));
			add2 = subjMat(s)-sum(vals);
			vals(1:add2) = [vals(1:add2)+1];
			siteDist{s} = vals;
			siteMat(s) = numel(vals);
		end
            end

            %% RUN ALGORITHM
            %%%The following code runs distributed ICA by performing the following
            % steps :
            % 1. Loading data into memory
            % 2. Masking Real Data
            % 3. Generating Simulated parts for Hybrid datasets and mixing with
            % real data
            % 4. Running Experiments for each desired subj/site combination

	    %%% 1. Loading Data
	    % For the code-release, it is assumed that data is input to the function
 
            %%% 2. Masking Real Data
            % for simulated 2-d data, just generate a circular mask
	    if flag_mask && load_mask == 0
                   mask =  compute_mask('maskType','mean','flag_transpose',0,datasets{:});
	    elseif load_mask ~= 0
		mask = load(load_mask);
		mask = mask.mask;
            end
            % the masking process is done by apply_mask
            % for more info on usage use 'help apply_mask'
            
            fprintf('\tMasking data outside of head...\n');
            if flag_mask
            	jump = compute_jump(length(datasets), 256, 2);	
		for i = 1:jump:length(datasets)
			m = i:i+(jump-1);
			fprintf('Masking %i through %i\n',m(1),m(end));
			parfor j = m
			datasets{j} = (apply_mask(datasets{j},'mask',mask','demean',0))';
		end
	    end
	end

	%% 4. Running experiments for each desired subj/site combination
	%%% this is the actual running of each set-up experiment
	% each experiment runs through the following steps for each subj/site
	% combination:
	% 4.1: set up run name and saving for the specific run
	% 4.2: distributing datasets across simulated nodes
	% 4.3: removing global mean across nodes
	% 4.4: performing local PCA 
	% 4.5: performing global PCA
	% 4.6: normalization
	% 4.7: distributed joint ICA
	% 4.8: Back-Reconstruction
	% 4.9: Zscoring
	% 4.10: Correlation with Groundtruth and Hungarian Algorithm
	% 4.11: ISI Computation
	% 4.12: refill mask
	% 4.13: Save results for run in .mat files
	% 4.14: Save spatial maps in .nii file
	for s = 1:numel(subjMat) %iterating over the subjMat seems ok
	%%% 4.1 Set the runname- for display, and saving
	% runname eg : s32-n2-nc52-r1
	runName = ['s' int2str(subjMat(s)) '-n' int2str(siteMat(s)) '-nc' int2str(nc) '-r' int2str(r)];

	%%% Varyon, Child and Nii directory (saved in the parent directory created
			%%% above) is just the runname, the nii for the spatial maps are
	%%% saved in the nii folder
	if ispc
	child_directory = sprintf('%s\%s',parent_directory,runName);
if iscell(varyS)
	child_directory = [child_directory '\' sprintf('%s-%s',varyOn,mat2str(varyS{v}))];
	else
	child_directory = [child_directory '\' sprintf('%s-%.3f',varyOn,varyS(v))];
	end
	else
	child_directory = sprintf('%s/%s',parent_directory,runName);
	if ~strcmp(varyOn,'STOP')
if iscell(varyS)
	child_directory = [child_directory '/' sprintf('%s-%s',varyOn,mat2str(varyS{v}))];
	else
	child_directory = [child_directory '/' sprintf('%s-%.3f',varyOn,varyS(v))];
	end
	end
	end
	if ~exist(child_directory,'dir')
	mkdir(child_directory);
	end


	%Checking to make sure neither is zeroed by accident
if (subjMat(s) ~= 0 && siteMat(s) ~= 0)    

	fprintf(['********************\nStarting  ' runName '\n********************\n\n']);
	fprintf('Performing Distributed Analysis...\n');

	%%% step 4.2: distributing sets across simulated nodes
	% datasets simulates sites- one element of the cell
	%represents a local dataset
	nodeDecomp = cell(1,siteMat(s)); %decomposed u_i

	% Sets are distributed across "nodes" according to
	% siteDistribution variable
	fprintf('Distributing data across sites...\n');
	% Instead of actually copying the distributed sets into cells,
	% we can save way more local memory by just saving the indices
	% which accord with the correct site distributions.
	dataIndex = zeros(2,siteMat(s));
	for i = 1:siteMat(s)
if ~iscell(siteDist)
	dataIndex(1,i) = (1 +(i-1)*siteDist(s));
	dataIndex(2,i) = i*siteDist(s);
	else
	sD = siteDist{s};
if i > 1
dataIndex(1,i) = (1+dataIndex(2,i-1));
dataIndex(2,i) = dataIndex(1,i) + sD(i) -1;
else
dataIndex(1,i) = 1;
dataIndex(2,i) = sD(1);
end

end
end

fprintf('Global set is size %d x %d\n',VOX(1),sum(TIME));

%%% step 4.3:  Remove mean... the global mean is computed using the method
% used by Bai, Chan, and Luk in their dPCA algorithm.
if flag_demean
fprintf('Demeaning...\n');
mu_n = 0;
mu_ni = 0;
mu = 0;
if meanDim == 2
	% Compute Local Means
for i = 1:siteMat(s)
	localMu = mean(cell2mat(datasets(dataIndex(1,i):dataIndex(2,i))),meanDim);
	mu_ni = size(cell2mat(datasets(dataIndex(1,i):dataIndex(2,i))),meanDim);
	mu = mu + localMu*mu_ni;
	mu_n = mu_n + mu_ni;
	end

	mu = mu/mu_n;
	% Demean Local Datasets according to Global Mean
	for i = 1:siteMat(s)
for j = 1:abs(dataIndex(2,i)-dataIndex(1,i))
	meanSub = repmat(mu,1,size(datasets{j},meanDim));
	datasets{j} = datasets{j} - meanSub;
	end
	end
	clearvars mu localMu
	else
for i = 1:siteMat(s)
	datasets{i} = datasets{i} - repmat(mean(datasets{i},meanDim),size(datasets{i},meanDim),1);
	end
	end
	end

	fprintf('Performing local pca...\n');
	%%% step 4.4: first level of PCA performed on each node in time
	%%% dimension
	Y = cell(1,siteMat(s));
	eigOpts.tol = 10^-6;
	eigOpts.disp = 0;
	eigOpts.maxit = 1000;
	jump = compute_jump(siteMat(s), 256, 2);
                      for i = 1:jump:length(Y)
                           m = i:(i+jump-1);
                           parfor j = m      
                               % Compute local PCA using Eigs function
				if flag_local_pca
                                	dd = cell2mat(datasets(dataIndex(1,j):dataIndex(2,j)));
                                	[u_i,~] = eigs(dd'*dd,k,'LM',eigOpts);
                                	% Project into local space
                                	Y{j} =  dd * u_i;
				else
					Y{j} = cell2mat(datasets(dataIndex(1,j):dataIndex(2,j)));
				end
                           end
                      end

                    if flag_verbose
                        fprintf('\n');
                    end

                    %%% step 4.5: distributed PCA - Rogers' method
		    if flag_dpca
                    	fprintf('\t\t\tPerforming distributed pca...\n');
                    	V = dpca(Y,nc);
                    	Y_ = cell(1,siteMat(s));
                    	for i = 1:siteMat(s)
                            Y_{i} = V'*cell2mat(datasets(dataIndex(1,i):dataIndex(2,i)));
                    	end
		    else
			Y_ = Y;
		    end
                    %%% step 4.6: Normalization
                    if  flag_normalize
                        fprintf('\t\t\tNormalizing data...\n');
                        % Normalize data globally (global data should havif mod(i,2) == 0, std dev = 1)
                        nm = diag(std([Y_{:}],[],2).^(-1)); % Normalizing factor for each component
                        parfor i = 1:siteMat(s)
                            Y_{i} = nm*Y_{i};   % Entire reduced, normalized, demeaned dataset
                            if flag_verbose
                                fprintf('/%i',i);
                            end
                        end
                    end


                    %%% step 4.7: distributed joint ICA 
		    if flag_djica
                    	fprintf('\t\t\tPerforming distributed ica...\n');
                    	% see help djica for more info
                   
                    	[W,sphere_n,activations_n,bias,lrates,changes] = djica(Y_{:});

                    	%%% step 4.8: Back- Reconstructing
                    	fprintf('\t\tbackreconstructing...\n');
                    	Shat = cell(1,siteMat(s));
                    	for i = 1:siteMat(s)            
                        	Shat{i} = W*nm*(V')*cell2mat(datasets(dataIndex(1,i):dataIndex(2,i)));
                    	end

                    	% Individual subjects need to be parsed out from the locally
                    	% stacked data
                    	Shat_ = cell(1,subjMat(s));
                    	if iscell(siteDist)
                        	sD = siteDist{s};
                    	else 
                        	sD = repmat(siteDist,1,subjMat(s));
                    	end
                    	for i = 1:siteMat(s)
                        	Shat_(dataIndex(1,i):dataIndex(2,i)) = mat2cell(Shat{i},nc,TIME(dataIndex(1,i):dataIndex(2,i)));
                    	end
                    	% Spatial Map is computed using the pseudo-inverse
                    	SMhat = pinv(W*(V'));
                    	% Save all of this at this point, in case of crash.
                    	save([child_directory '/ica_materials.mat'],'W','nm','V');

                  	  %%% Step 4.9: zscore maps
			if flag_zscore
                    		fprintf('\t\tZscoring maps...\n');
                    		for i = 1:nc
                        		SMhat(i,:) = zscore(SMhat(i,:));
                    		end
			end

                   	 %%% Step 4.10: GT Correlation and Hungarian Algorithm
                    	% only correlate with the ground-truth if we have simulated
                    	% components with which to correlate
                    
                    	%%% step 4.11 ISI computation
			ISI = 0;
                    	if flag_isi	
				try
					load('ground_truth.mat');
       		               		 subject_corr = cell(1,length(Shat_));
                        		for i = 1:length(Shat_)
                            			Shat_{i} = Shat_{i}'
                            			subject_corr{i} = abs(corr(Shat_{i},times{i}'))';
                        		end
                       			A = apply_mask(sims{1},'mask',mask','demean',0)';
                        		H = hungarian_algorithm(SMhat',subject_corr{i});
                      			P = H*A';
                     			 % this uses compute_isi.m, see 'help cocompute_isi' for more info
                        		ISI  = compute_isi(P,nc);
                        		fprintf('\t\tSpatial Map ISI = %f\n',ISI(1,s));

                 		       %%% Saving
                        		save([child_directory '/isi.mat'],'ISI');
                        		save([child_directory '/subject_corr.mat'],'subject_corr');
				catch E
					fprintf('EXCEPTION %s\nFAILED TO COMPUTE ISI. Perhaps you are using custom data, or ground_truth.at is not available?\n', E.message);
				end
			end
		    end
                    %%% step 4.12 : filling mask
                    if flag_mask 
                        save([child_directory '/mask.mat'],'mask');
                    end

                    %%% step 4.13 : saving results in .mat files
		    if flag_djica
                    	save([child_directory '/convergence.mat'],'changes','lrates'); % save timings
                   	save([child_directory '/SMhat.mat'],'SMhat'); % save the spatial map at this point...
                    	save([child_directory '/IC.mat'],'Shat_'); % save subject-specific timecourses
		    end
                end

                fprintf(['********************\nDone with ' runName '\n********************\n\n']);  
            end
       end
    end
    %% Output Post-Processing
    fprintf('DONE WITH ALL TESTS!!!!\n');
end

%djica() - Perform Independent Component Analysis (ICA) decomposition
%            of psychophysiological data using the infomax ICA algorithm of 
%            Bell & Sejnowski (1995) with the natural gradient feature 
%            of Amari, Cichocki & Yang in a distributed framework. 
%
% Usage:
%      simply >> [weights] = djica(data1,data2,data3,...,dataN);
%       or    
%        else >> [weights,sphere,bias,lrates] ...
%                                 = djica(data1,data2,data3,...,dataN,'Key1',Value1',...);
%
% Input_Variable:
%
% datai  - the ith (channels x frames) local datasets 
%
% Optional_Keywords       Keyword_Values                  Default_Values
%
% 'ncomps'    = [N] number of ICA components to compute (default -> chans)
%               using rectangular ICA decomposition
% 'sphering'  = ['on'/'off'] flag sphering of data      (default -> 'on')
% 'weights'   = [W] initial weight matrix               (default -> eye())
%                            (Note: if 'sphering' 'off', default -> spher())
% 'lrate'     = [rate] initial ICA learning rate (<< 1) (default -> heuristic)
% 'block'     = [N] ICA block size (<< datalength)      (default -> heuristic)
% 'anneal'    = annealing constant (0,1] (defaults -> 0.90, or 0.98, extended)
%                         controls speed of convergence
% 'annealdeg' = [N] degrees weight change for annealing (default -> 70)
% 'stop'      = [f] stop training when weight-change < this (default -> 1e-6)
% 'maxsteps'  = [N] max number of ICA training steps    (default -> 512)
% 'bias'      = ['on'/'off'] perform bias adjustment    (default -> 'on')
% 'posact'    = make all component activations net-positive(default 'on'}
% 'verbose'   = give ascii messages ('on'/'off')        (default -> 'on')
%
% Output_Variables [RO = output in reverse order of projected mean variance 
%                        unless starting weight matrix passed ('weights' above)]
%
% weights     = Global ICA weight matrix (comps,chans)     [RO]
% sphere      = cell containing the data sphering matrices (chans,chans) = spher(data)
%               Note: unmixing_matrix = weights*sphere {sphering off -> eye(chans)}
% bias        = vector of final (ncomps) online bias [RO]    (default = zeros())
% lrates      = vector of learning rates used at each training step
%
%
%
% For more information:
% http://www.cnl.salk.edu/~scott/icafaq.html - FAQ on ICA/EEG
% http://www.cnl.salk.edu/~scott/icabib.html - mss. on ICA & biosignals
% http://www.cnl.salk.edu/~tony/ica.html - math. mss. on ICA, with kernal code
%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Edit history %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	Based on:
%  runica()  - by Scott Makeig with contributions from Tony Bell, Te-Won Lee 
%              Tzyy-Ping Jung, Sigurd Enghoff, Michael Zibulevsky et al.
%                            CNL / Salk Institute 1996-99
%
%		Built based on icatb_runica - extraneous functions removed and cleaned up
%		Distributed ICA for two nodes
%		Distributed ICA for n nodes
%		Preprocessing move to outside function
%		Parallel processing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

function [weights,sphere_n,activations_n,bias,lrates,changes,runtimes_n] = djica(varargin)

if length(varargin) < 1
   help djica
   return
end

%argument processing function included below
[keys,vals,datasets] = process_args(varargin,1);
%% Other preprocessing

% chans_n{i} are the rows (components) of the local dataset i, and frames_n{i} is the
% number of observations
[chans_n frames_n] = cellfun(@size,datasets,'uniformOutput',1);
runtimes_n = zeros(1,length(datasets));
datalength_n = frames_n;
if arrayfun(@(x) x<2,chans_n) | arrayfun(@(x,y) x<y,frames_n,chans_n) 
    fprintf('\ndjica() - data size too small.\n\n');
    fprintf('\ndjica() - flipping matrix...\n\n');
    datasets = cellfun(@transpose,datasets,'uniformoutput',false);
end

if range(chans_n) ~= 0 
	fprintf('\ndjica() - differing number of components.\n\n');
	return
end

%
%% Declare defaults used below %%%%%%%%%%%%%%%%%%%%%%%%
%
MAX_WEIGHT           = 1e8;       % guess that weights larger than this have blown up
DEFAULT_STOP         = 0.000001;  % stop training if weight changes below this
DEFAULT_ANNEALDEG    = 60;        % when angle change reaches this value,
DEFAULT_ANNEALSTEP   = 0.90;      %     anneal by multiplying lrate by this

DEFAULT_MAXSTEPS     = 1024;       % stop training after this many steps 

DEFAULT_BLOWUP       = 1000000000.0;   % = learning rate has 'blown up'
DEFAULT_BLOWUP_FAC   = 0.8;       % when lrate 'blows up,' anneal by this fac
DEFAULT_RESTART_FAC  = 0.9;       % if weights blowup, restart with lrate
% lower by this factor
MIN_LRATE            = 0.000001;  % if weight blowups make lrate < this, quit
MAX_LRATE            = 0.1;       % guard against uselessly high learning rate
DEFAULT_LRATE        = (log(min(chans_n)).^(-1)).*0.015; % per-site learning rate 
% heuristic default - may need adjustment
%   for large or tiny data sets!
DEFAULT_BLOCK        = floor((min(frames_n)./20).^(1/2));  % heuristic default - uses lowest frame-value
% - may need adjustment!

DEFAULT_SPHEREFLAG   = 'on';      % use the sphere matrix as the default
%   starting weight matrix
DEFAULT_POSACTFLAG   = 'on';      % use posact()
DEFAULT_VERBOSE      = 1;         % write ascii info to calling screen
DEFAULT_BIASFLAG     = 1;         % default to using bias in the ICA update rule
DEFAULT_OPTIMIZATION = 'stochastic';
%                                 
%%%%%%%%%%%%%%%%%%%%%%% Set up keyword default values %%%%%%%%%%%%%%%%%%%%%%%%%
%
epochs = 1;							 % do not care how many epochs in data

sphering   = DEFAULT_SPHEREFLAG;     % default flags
posactflag = DEFAULT_POSACTFLAG;
verbose    = DEFAULT_VERBOSE;

block      = DEFAULT_BLOCK;          % heuristic default - may need adjustment!
lrate      = DEFAULT_LRATE;
annealdeg  = DEFAULT_ANNEALDEG;
annealstep = 0;                      % defaults declared below
nochange   = DEFAULT_STOP;
maxsteps   = DEFAULT_MAXSTEPS;

weights    = 0;                      % defaults defined below

ncomps     = chans_n(1);
biasflag   = DEFAULT_BIASFLAG;

wts_blowup = 0;                      % flag =1 when weights too large
wts_passed = 0;                      % flag weights passed as argument
optimization = DEFAULT_OPTIMIZATION;
%
%% Collect keywords and values from argument list %%%%%%%%%%%%%%%
%

for i = 1:length(keys)% for each Keyword
   Keyword = keys{i};
   Value = vals{i};
   if ~ischar(Keyword)
      fprintf('djica(): keywords must be strings')
      return
   end
   Keyword = lower(Keyword); % convert upper or mixed case to lower
   
   if strcmp(Keyword,'weights') | strcmp(Keyword,'weight')
      if ischar(Value)
         fprintf(...
            'djica(): weights value must be a weight matrix or sphere')
         return
      else
         weights = Value;
         wts_passed =1;
      end
   elseif strcmp(Keyword,'posact') 
      if ~ischar(Value)
         fprintf('djica(): posact value must be on or off')
         return
      else 
         Value = lower(Value);
         if ~strcmp(Value,'on') & ~strcmp(Value,'off'),
            fprintf('djica(): posact value must be on or off')
            return
         end
         posactflag = Value;
      end
   elseif strcmp(Keyword,'lrate')
      if ischar(Value)
         fprintf('djica(): lrate value must be a number')
         return
      end
      lrate = Value;
      if lrate>MAX_LRATE | lrate <0,
         fprintf('djica(): lrate value is out of bounds'); 
         return
      end
      if ~lrate,
         lrate = DEFAULT_LRATE;
      end
   elseif strcmp(Keyword,'block') | strcmp(Keyword,'blocksize')
      if ischar(Value)
         fprintf('djica(): block size value must be a number')
         return
      end
      block = Value;
      if numel(block) == 1
          block = repmat(block,1,length(datasets));
      end
      if ~block,
         block = DEFAULT_BLOCK; 
      end
   elseif strcmp(Keyword,'stop') | strcmp(Keyword,'nochange') ...
         | strcmp(Keyword,'stopping')
      if ischar(Value)
         fprintf('djica(): stop wchange value must be a number')
         return
      end
      nochange = Value;
   elseif strcmp(Keyword,'maxsteps') | strcmp(Keyword,'steps')
      if ischar(Value)
         fprintf('djica(): maxsteps value must be an integer')
         return
      end
      maxsteps = Value;
      if ~maxsteps,
         maxsteps   = DEFAULT_MAXSTEPS;
      end
      if maxsteps < 0
         fprintf('djica(): maxsteps value must be a positive integer')
         return
      end
   elseif strcmp(Keyword,'anneal') | strcmp(Keyword,'annealstep')
      if ischar(Value)
         fprintf('djica(): anneal step constant must be a number (0,1)')
         return
      end
      annealstep = Value;
      if annealstep <=0 | annealstep > 1,
         fprintf('djica(): anneal step value must be (0,1]')
         return
      end
   elseif strcmp(Keyword,'annealdeg') | strcmp(Keyword,'degrees')
      if ischar(Value)
         fprintf('djica(): annealdeg value must be a number')
         return
      end
      annealdeg = Value;
      if ~annealdeg,
         annealdeg = DEFAULT_ANNEALDEG;
      elseif annealdeg > 180 | annealdeg < 0
         fprintf('djica(): annealdeg value is out of bounds [0,180]')
         return
         
      end
   elseif strcmp(Keyword,'sphering') | strcmp(Keyword,'sphereing') ...
         | strcmp(Keyword,'sphere')
      if ~ischar(Value)
         fprintf('djica(): sphering value must be on, off, or none')
         return
      else 
         Value = lower(Value);
         if ~strcmp(Value,'on') & ~strcmp(Value,'off') & ~strcmp(Value,'none'),
            fprintf('djica(): sphering value must be on or off')
            return
         end
         sphering = Value;
      end
   elseif strcmp(Keyword,'bias')
      if ~ischar(Value)
         fprintf('djica(): bias value must be on or off')
         return
      else 
         Value = lower(Value);
         if strcmp(Value,'on') 
            biasflag = 1;
         elseif strcmp(Value,'off'),
            biasflag = 0;
         else
            fprintf('djica(): bias value must be on or off')
            return
         end
      end
   elseif strcmp(Keyword,'verbose') 
      if strcmp(Value,'on') || Value == 1,
         verbose = 1; 
      elseif strcmp(Value,'off') || Value == 0,
         verbose = 0; 
      else
         fprintf('djica(): verbose flag value must be on or off')
         return
      end
   else
      fprintf('djica(): unknown flag')
      return
   end
end
%
%% Initialize weights, etc. %%%%%%%%%%%%%%%%%%%%%%%%
%
if ~annealstep,
      annealstep = DEFAULT_ANNEALSTEP;     % defaults defined above
end % else use annealstep from commandline

if ~annealdeg, 
   annealdeg  = DEFAULT_ANNEALDEG - momentum*90; % heuristic
   if annealdeg < 0,
      annealdeg = 0;
   end
end
if ncomps >  chans_n(1) | ncomps < 1
   fprintf('djica(): number of components must be 1 to %d.\n',chans);
   return
end

%
%% Check keyword values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if sum(frames_n)<sum(chans_n),
   fprintf('djica(): data length %d < data channels %f!\n',frames,chans)
   return
elseif any(block < 2),
   fprintf('djica(): block size %d too small!\n',block)
   return
elseif any(block > min(frames_n)), 
   fprintf('djica(): block size exceeds data length!\n');
   return
elseif floor(epochs) ~= epochs,
   fprintf('djica(): data length is not a multiple of the epoch length!\n');
   return
end;
%
%% Process the data %%%%%%%%%%%%%%%%%%%%%%%%%%
%
if verbose,
	for i = 1:length(datasets)
		fprintf( ...
		   '\nInput data size for set %d [%d,%d] = %d channels, %d frames.\n', ...
		   i,chans_n(i),frames_n(i),chans_n(i),frames_n(i));
       fprintf('Initial learning rate for set %d will be %g, block size %d.\n',i,lrate,block);
	end
   
   fprintf( ...
      'Learning rate will be multiplied by %g whenever angledelta >= %g deg.\n', ...
      annealstep,annealdeg);
   fprintf('Training will end when wchange < %g or after %d steps.\n', ...
      nochange,maxsteps);
   if biasflag,
      fprintf('Online bias adjustment will be used.\n');
   else
      fprintf('Online bias adjustment will not be used.\n');
   end
end
%

%% Perform Local Sphering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
sphere_n = {};

if strcmp(sphering,'on'), %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if verbose,
      fprintf('Computing the sphering matrix...\n');
   end
   sphere_n = cellfun(@(x) 2.0*inv(sqrtm(cov(x'))),datasets,'uniformoutput',false);
   datasets = cellfun(@(x,y) x*y,sphere_n,datasets,'uniformoutput',false);
   if ~weights,
		if verbose,
		   fprintf('Starting weights are the identity matrix ...\n');
		end
   	weights = eye(ncomps,chans_n(1)); % begin with the identity matrix
   else % weights given on commandline
      if verbose,
         fprintf('Using starting weights named on commandline ...\n');
      end
   end
   if verbose,
      fprintf('Sphering the data ...\n');
   end   
elseif strcmp(sphering,'off') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if ~weights
      if verbose,
         fprintf('Using the sphering matrix as the starting weight matrix ...\n');
         fprintf('Returning the identity matrix in variable "sphere" ...\n');
      end
      sphere_n = cellfun(@(x) 2.0*inv(sqrtm(cov(x'))),datasets,'uniformoutput',false);
      weights = sphere_n{1};
      for i = 1:length(sphere_n)
          sphere_n{i} = eye(chans_n(i));
      end
   else % weights ~= 0
		if verbose,
 	     fprintf('Using starting weights named on commandline ...\n');
	     fprintf('Returning the identity matrix in variable "sphere" ...\n');
        end
        for i = 1:length(sphere_n)
          sphere_n{i} = eye(chans_n(i));
        end
   end
elseif strcmp(sphering,'none')
      sphere_n = cell(1,length(datasets));
      for i = 1:length(datasets)
          sphere_n{i} = eye([chans_n(i) chans_n(i)]);
      end
     if ~weights
      if verbose,
         fprintf('Starting weights are the identity matrix ...\n');
         fprintf('Returning the identity matrix in variable "sphere_n" ...\n');
      end
%       weights = cellfun(@(x) eye(ncomps,x),chans_n,'uniformoutput',false);
      weights = eye([chans_n(1) chans_n(1)]);
      
   else % weights ~= 0
      if verbose,
         fprintf('Using starting weights named on commandline ...\n');
         fprintf('Returning the identity matrix in variable "sphere" ...\n');
      end
   end
   if verbose,
      fprintf('Returned variable "sphere" will be the identity matrix.\n');
   end
end
%
%% Initialize ICA training %%%%%%%%%%%%%%%%%%%%%%%%%
%

delta=zeros(1,chans_n(1)*ncomps);
changes = [];
degconst = 180./pi;
startweights = weights;
startgrad = eye(ncomps,ncomps);
prevweights = startweights;
prevgrads = startgrad;
oldweights = startweights;
prevwtchange = zeros(chans_n(1),ncomps);
oldwtchange = zeros(chans_n(1),ncomps);
lrates = zeros(1,maxsteps);

%onesrow = ones(1,min(frames_n));
bias = zeros(ncomps,1);

bias_n = cell(1, length(datasets));
for i = 1:length(datasets)
    bias_n{i} = zeros(ncomps,1);
end

%
%% ICA training loop using the logistic sigmoid %%%%%%%%%%%%%%%%%%%
%
if verbose,
   fprintf('Beginning ICA training ...\n');
end
step=0;
laststep=0; 
blockno = 1;  % running block counter for kurtosis interrupts

permute_n = cell(1,length(datasets));
while step < maxsteps, %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for i = 1 :length(datasets)
           tic
    	   permute_n{i}=randperm(datalength_n(i)); % shuffle data order at each step, for each dataset
           runtimes_n(i) = runtimes_n(i) + toc;
   end
   %%%%%%%%% ICA Training Block %%%%%%%%%%%%%%%%%%%
		gradSum = zeros(chans_n(1),chans_n(1));
        biasSum = zeros(ncomps,1);
        grads = cell(1,length(datasets));
        biases = cell(1,length(datasets));
        for i = 1:length(datasets)
            grads{i} = gradSum;
            biases{i} = biasSum;
        end
        if strcmp(optimization, 'stochastic')
            lastt=fix((min(frames_n)/block-1)*block+1); %uses same heuristic as for determining block-size - minimum frame count
            onesrow = ones(1,block);
            BI=block*eye(ncomps,ncomps);
            for t=1:block:lastt,
                gradSum = zeros(chans_n(1),chans_n(1));
                biasSum = zeros(ncomps,1);
                grads = cell(1,length(datasets));
                biases = cell(1,length(datasets));
                for i = 1:length(datasets)
                    grads{i} = gradSum;
                    biases{i} = biasSum;
                end
                parfor i = 1:length(datasets)
                    if biasflag           
                      dat = datasets{i};
                      u=weights*dat(:,permute_n{i}(t:t+block-1)) + bias*onesrow;
                    else                                                             
                      u=weights*dat(:,permute_n{i}(t:t+block-1));
                    end                                                              

                      %%%%%%%%%%%%%%%%%%% Logistic ICA weight update %%%%%%%%%%%%%%%%%%%
                      y=1./(1+exp(-u));                                                %
                      grad=  lrate*(BI+(1-2*y)*u')*weights;                   %        
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                   if biasflag 
                         %%%%%%%%%%%%%%%%%%%%%%%% Logistic ICA bias %%%%%%%%%%%%%%%%%%%%%%%
                         bias_grad = lrate*sum((1-2*y)')'; % for logistic nonlin. %      
                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                              
                   end
                   grads{i} = grads{i} + grad;
                   biases{i} = biases{i} + bias_grad;
                end
                for i = 1:length(datasets)
                    gradSum = gradSum + grads{i};
                    biasSum = biasSum + biases{i};
                end
                weights = weights + gradSum; % weights is a global variable
                bias =  bias + biasSum; % bias is a global variable
            end
        else
            parfor i = 1:length(datasets) % for each local dataset
                tic
                onesrow = ones(1,block);
                BI=block*eye(ncomps,ncomps);
                lastt=fix((frames_n(i)/block-1)*block+1); %uses same heuristic as for determining block-size - minimum frame count
                for t=1:block:lastt,
                    if biasflag           
                      dat = datasets{i};
                      u=weights*dat(:,permute_n{i}(t:t+block-1)) + bias*onesrow;
                    else                                                             
                      u=weights*dat(:,permute_n{i}(t:t+block-1));
                    end                                                              

                      %%%%%%%%%%%%%%%%%%% Logistic ICA weight update %%%%%%%%%%%%%%%%%%%
                      y=1./(1+exp(-u));                                                %
                      grad=  lrate*(BI+(1-2*y)*u')*weights;                   %        
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                   if biasflag 
                         %%%%%%%%%%%%%%%%%%%%%%%% Logistic ICA bias %%%%%%%%%%%%%%%%%%%%%%%
                         bias_grad = lrate*sum((1-2*y)')'; % for logistic nonlin. %      
                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                              
                   end
                   grads{i} = grads{i} + grad;
                   biases{i} = biases{i} + bias_grad;
                end
                runtimes_n(i) = runtimes_n(i) + toc;
            end
            for i = 1:length(datasets)
              gradSum = gradSum + grads{i};
              biasSum = biasSum + biases{i};
            end
        end
      %%%%%%%%%%%%%%%%%%% Aggregation %%posactWeights%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if ~strcmp(optimization, 'stochastic')
          weights = weights + gradSum; % weights is a global variable
          bias =  bias + biasSum; % bias is a global variable
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      if max(max(abs(weights)))> MAX_WEIGHT
         wts_blowup = 1;
         change = nochange;
      end
      
   
   if ~wts_blowup
      oldwtchange = weights-oldweights;
      step=step+1; 
      %
      %%%%%%% Compute and print weight and update angle changes %%%%%%%%%
      %
      lrates(1,step) = mean(lrate);
      angledelta=0.;
      delta=reshape(oldwtchange,1,chans_n(1)*ncomps);
      change=delta*delta'; 
   end
   %
   %%%%%%%%%%%%%%%%%%%%%% Restart if weights blow up %%%%%%%%%%%%%%%%%%%%
   %
   if wts_blowup | isnan(change)|isinf(change),  % if weights blow up,
      fprintf('');
      step = 0;                          % start again
      change = nochange;
      wts_blowup = 0;                    % re-initialize variables
      blockno = 1;
      lrate = lrate*DEFAULT_RESTART_FAC; % with lower learning rate
      weights = startweights;            % and original weight matrix
	  grad = startgrad;
      gradSum = zeros(chans_n(1),chans_n(1));
      biasSum = zeros(ncomps,1);
      oldweights = startweights;            
      change = nochange;
      oldwtchange = zeros(chans_n(1),ncomps);
      delta=zeros(1,chans_n(1)*ncomps);
      olddelta = delta;
 
      prevweights = startweights;
      prevwtchange = zeros(chans_n(1),ncomps);
      lrates = zeros(1,maxsteps);
      bias = zeros(ncomps,1);
		for i = 1:length(datasets)
			bias_n{i} = zeros(ncomps,1);
		end

      if any(lrate> MIN_LRATE)
         r = rank(datasets{1});
         if r<ncomps
            fprintf('Data has rank %d. Cannot compute %d components.\n',...
               r,ncomps);
            return
         else
            fprintf(...
               'Lowering learning rate to %g and starting again.\n',lrate);
         end
      else
         fprintf( ...
            'runica(): QUITTING - weight matrix may not be invertible!\n');
         return;
      end
   else % if weights in bounds 
      %
      %%%%%%%%%%%%% Print weight update information %%%%%%%%%%%%%%%%%%%%%%
      %
      if step> 2 
         angledelta=acos((delta*olddelta')/sqrt(change*oldchange));
      end
      if verbose,
         if step > 2, 
           
            fprintf(...
               'step %d - lrate %5f, wchange %7.6f, angledelta %4.1f deg\n', ...
               step,lrate,change,degconst*angledelta);
            
         
         else
            fprintf(...
               'step %d - lrate %5f, wchange %7.6f\n',step,lrate,change);
         end % step > 2
      end; % if verbose
      %
      %%%%%%%%%%%%%%%%%%%% Save current values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %
      changes = [changes change];
      oldweights = weights;
      %
      %%%%%%%%%%%%%%%%%%%% Anneal learning rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %
      if degconst*angledelta > annealdeg,  
         lrate = lrate*annealstep;          % anneal learning rate
         olddelta   = delta;                % accumulate angledelta until
         oldchange  = change;               %  annealdeg is reached
      elseif step == 1                     % on first step only
         olddelta   = delta;                % initialize 
         oldchange  = change;               
      end
      %
      %%%%%%%%%%%%%%%%%%%% Apply stopping rule %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %
      if step >2 & change < nochange,      % apply stopping rule
         laststep=step;            
         step=maxsteps;                  % stop when weights stabilize
      elseif change > DEFAULT_BLOWUP,      % if weights blow up,
         lrate=lrate*DEFAULT_BLOWUP_FAC;    % keep trying 
      end;                                 % with a smaller learning rate
   end; % end if weights in bounds
   
end; % end training %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~laststep
   laststep = step;
end;
lrates = lrates(1,1:laststep);           % truncate lrate history vector
%
%% Orient components towards positive activation %%%%%%%%%%%
%
if strcmp(posactflag,'on')
    activations_n = cell(1,length(datasets));
	for i = 1:length(datasets)
        [activations_n{i}, winvout, weight_p] = icatb_posact(datasets{i}, weights);
		  
    end
   % changes signs of activations and weights to make activations
   % net rms-positive
else
    activations_n = cell(1,length(datasets));
    for i = 1:length(datasets)
        activations_n{i} = weights*datasets{i};	  
    end
end
%

%
%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
return

if nargout > 6
   u=weights*data + bias*ones(1,frames);      
   y = zeros(size(u));
   for c=1:chans
      for f=1:frames
         y(c,f) = 1/(1+exp(-u(c,f)));
      end
   end
end
end

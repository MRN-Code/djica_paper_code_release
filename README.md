# djICA paper code release

This repository contains the MATLAB code used in the 2018 paper "Decentralized Temporal Independent Component Analysis: Leveraging fMRI Data in Collaborative Settings". 

## Running the Pipeline

We have included some sample simulated data with 512 simulated subjects. The pipeline can be run for this data set by running the following

  ```MATLAB
    load datasets
    ISI = djica_pipeline(datasets);
  ```
  
  which runs 512 subjects simulated over 64 sites. 

The pipeline saves results in folders of the format (e.g. for 512 subjects 64 sites 20 ICs):

  ./results/s512-n64-nc20-r1/
  
The following objects are saved in a full run:

  * 

The code is written to run with MATLAB's distributed computing toolbox, but will run in serial mode if the toolbox is
not loaded beforehand. Running with the parallel toolbox decreases runtime of all stages significantly.

To run specific parts of the pipeline, you can pass flags. For example,

  ```MATLAB
  djica_pipepline(datasets, 'flag_local_pca', 1, 'flag_dpca', 0, 'flag_djica', 0);
  ```
runs only the local PCA stage on the input data.

### Using custom data

To run on your own data set, assemble your data into a cell, where each entry of the cell
is an individual subject. If you skip the local PCA step, the pipeline will assume that each
entry has already been processed with local PCA, and each entry represents a site.

This cell can then be passed into the djica_pipeline.

If you want to compute ISI with Real Data, or custom simulated data, you will need to save your ground truth,
or pseudo-ground-truth into ground_truth.mat. There are two variables stored here, 'sims', which
stores ground truth ICxVOXEL spatial maps, and 'times' which stores ground truth TIMExIC timecourses.

For the real data case, you will need to back-reconstruct these yourself. If the proper ground-truth is not
provided, the pipeline will still run, but ISI will not be computed.

### Additional Options and Flags

Other than the flags mentioned above, the pipeline has a large number of keywords and variables that can be controlled from the MATLAB pipeline. Type 'help djica_pipeline', orlook at the first section of djica_pipeline.m for more detailed arguments that can be passed to the function.

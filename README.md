# SIDELINE
The SIngle-cell, DEep-Learning gene regulatory network INferencE (SIDELINE) method is a computational method designed to predict gene regulatory networks (GRNs) using paired datasets of case versus control experiments. After constructing a gene co-differential expression network, SIDELINE employs cell-based pseudotiming, an attention-based convolutional neural network method, and permutation-based significance testing for inferring GRNs among gene modules. 

Corresponding Paper: ADD PAPER LINK HERE 
![MyImage](SIDELINE.png)

## Functionality
* Bullet main uses for SIDELINE 


## Prerequisites 
### General 
- Python >= 3.5
- [PyTorch](https://pytorch.org/get-started/locally/)
- Optional: CUDA 
### Required Python Packages 
- numpy
- pandas
- random
- heapq
- copy
- os
- sys
- matplotlib
- pylab
- networkx
- argparse
- bambi
- scipy
- scanpy
- time
- arviz
- math
- argparse
- warnings 
- shutil
### Data
Required: Two scRNA-seq datasets, one case and one control
Format: A non-normalized CSV file with genes as rows and cells as columns. **Genes must be in first column of CSV file** 

*Note: Case and control datasets should either be the same cell type and two different experimental conditions OR the same experimental condition and two different cell types* 
##### Sample Datasets 
There are multiple sample datasets under the [Data folder](https://github.com/chenyongrowan/SIDELINE/tree/main/Data). 
1. The [ProstateCancer](https://github.com/chenyongrowan/SIDELINE/tree/main/Data/ProstateCancer) folder contains datasets for 4 patients. The cells were processed to now contain only endothelial cells. 
2. The [RemoteMemoryFormation](https://github.com/chenyongrowan/SIDELINE/tree/main/Data/RemoteMemoryFormation) folder contains preprocessed datasets from 2 papers. Both datasets contain only neurons. (FC = Fear Conditioned). The [Rao-Ruiz](https://github.com/chenyongrowan/SIDELINE/tree/main/Data/RemoteMemoryFormation/Rao-Ruiz) dataset is smaller, demonstrating functionality with a small cell count. 
### Running 
SIDELINE is set up as a single-line command with the following flags: 
**Required flags**
- -goi/--geneOfInterest:  One or more genes of interest. Separate multiple genes with a "+" (ex. Arc+Bdnf)
- -ctrl/--control:        Path to the csv file containing control cells. Provide file with cells as columns and genes as rows. The gene names should be the first column in the file. *The file must contain at least 10 cells*
- -exp/--experimental:    Path to the csv file containing the case cells. Provide file with cells as columns and genes as rows. The gene names should be the first column in the file. *The file must contain at least 10 cells*

**Optional flags**
- -p/--permutations:  Number of permutations to run. Default 100
- -top/--numTopGenes: Number of top correlated genes selected. Default 100
- -zero/--zeroThresh: Threshold for number of 0's tolerated for a gene. Default 0.30
- -s/--start:         Starting point for SIDELINE. Default 1 (Run SIDELINE and background). 2 runs only the background
- --cuda:             CUDA use on when flag included. Leave flag out if using CPU based discovery
- -o/--output:        Output directory name. Default is 'SIDELINE_Output' 

Give command to run from SIDELINE download for sample datasets


## Paper
Corresponding Paper (peer-reviewed, open access): ADD PAPER NAME AND LINK 

Please cite this paper when using SIDELINE. 


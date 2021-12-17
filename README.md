[![DOI:10.5281/zenodo.5786211](https://zenodo.org/badge/DOI/10.5281/zenodo.5786211.svg)](https://doi.org/10.5281/zenodo.5786211)

Introduction
------------------------

README file for codes to reproduce the three downstream analyses (such as *differential expression analysis*, *cell clustering analysis* and *pseudotime analysis*) in the paper **"Imputing dropouts for single-cell RNA sequencing based on multi-objective optimization"**. The method `scMOO` developed in the paper can be found at: [https://github.com/Zhangxf-ccnu/scMOO](https://github.com/Zhangxf-ccnu/scMOO).   
To reproduce the *masking* and *down-sampling* experiments, please refer to Chen and Zhou (2018): [Vpaper2018](https://github.com/ChenMengjie/Vpaper2018).


Contents of this archive
------------------------
This archive contains:   
 
(1) **Datasets**: subdirectory that contains four preprocessed datasets: bulk data and single-cell data of **H1\_DEC** (such as **H1\_DEC\_bulk** and **H1\_DEC\_sc**), **PBMC\_CL** and **Deng** datasets, which can be used to reproduce the three downstream analyses (such as *differential expression analysis*, *cell clustering analysis* and *pseudotime analysis*) respectively.   
Note that **H1\_DEC** and **PBMC\_CL** datasets have been preprocessed using Seurat v3.2 to contain 2,000 highly variable genes, and **Deng** dataset has been preprocessed by filtering out genes expressed in less than 10% of cells.


(2) **DE\_analysis**: subdirectory that contains three R codes to reproduce the differential expression analysis.    

`Step1_Prepare_datasets.R`: After downloading the bulk data and single-cell data of **Cell Type (GSE75748)** from GEO website, selecting 6 pairs of cell subpopulations including DEC, then using Seurat v3.2 to select 2,000 highly variable genes of both the bulk data and single-cell data.    

`Step2_edgeR.R`: Using `edgeR` to identify differential expression genes (DEGs) between pairs of cell subpopulations.  

`Step3_Compute_AUCscores_SpearmanCorrelation.R`: Computing AUC scores and Spearman correlation coefficient.

(3) **Cell\_clustering**: subdirectory that contains two R codes to reproduce the cell clustering analysis. 

`Step1_Select_HVGs.R`: Using Seurat v3.2 to select 2,000 highly variable genes of the single-cell data.    

`Step2_SC3_Seurat.R`: Using `SC3` and `Seurat` to carry out cell clustering analysis. And `ARI` and `NMI` are used to evaluate the consistency between the results of `SC3` or `Seurat` and the reference labels of cells.   

(4) **Cell\_trajectory\_inference**: subdirectory that contains three R codes to reproduce the pseudotime analysis. 

`Step1_Preprocess.R`: Obtaining the preporcessed **Deng** dataset with settings `percent=0.1` and `preprocess.only=TRUE`. 

`Step2_Monocle2.R`: Using `Monocle2` to carry out pseudotime analysis. The function returns `POS` and `Kendall's rank correlation` scores.   

`Step3_Trajectory_inference.R`: Running `Monocle2` on preprocessed dataset with corresponding setting `cellLabels`.


Tutorial
------------------------
A tutorial with example of cell clustering analysis at the same time illustrating the usage of `scMOO` is available at:
[scMOO-tutorial](https://github.com/Zhangxf-ccnu/scMOO) 


Contact
------------------------
Please do not hesitate to contact Miss **Ke Jin** [kej13@mails.ccnu.edu.cn](kej13@mails.ccnu.edu.cn) or Dr. **Xiao-Fei Zhang** [zhangxf@mail.ccnu.edu.cn](zhangxf@mail.ccnu.edu.cn) to seek any clarifications regarding any contents or operation of the archive.

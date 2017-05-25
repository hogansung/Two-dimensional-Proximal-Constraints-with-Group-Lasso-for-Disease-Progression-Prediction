# Two-dimensional Proximal Constraints with Group Lasso for Disease Progression Prediction

## Methodology
In this paper, I mainly contribute in extending multitask learning models with one-dimensional constraint into model with two-dimensional ones.
- Extension From 1D-TGL to 2D-TGL and 2D-TGL+
- Extension From 1D-cFSGL to 2D-cFSGL and 2D-cFSGL+


## Dataset
In this experiment, I use dataset provided by Alzheimer’s Disease Neuroimaging Initiative (ADNI), whose information can be found in their [website](http://adni.loni.usc.edu/).


### How To Obtain Dataset?
Since the dataset is preserved for request-only. unfortunately, I cannot provide the dataset publicly. To apply for the dataset, please refer to this [link](http://adni.loni.usc.edu/data-samples/access-data/).


### What Is Included In The Dataset?
To fully reproduce my experiment results, a list of dataset is required. Files provide two kinds of features as follows.
1. MRI features
    - UCSFFSL_02_01_16.csv
2. META features 
    - ADASSCORES.csv
    - CDR.csv
    - FAQ.csv
    - GDSCALE.csv
    - LABDATA.csv
    - MODHACH
    - NEUROBAT.csv
    - PTDEMOG.csv
    - MMSE.csv

A detailed information for META features are shown as below.
![](img/META_feature_info.png)


### How To Preprocess Dataset?
Once all the dataset is downloaded, please make sure to put all those downloaded files under the folder ['dat/Origin'](dat/Origin) and create folder ['dat/TGL_sgl/Longitudinal/'](dat/TGL_sgl/Longitudinal/) and ['dat/TGL_mtl/Longitudinal/'](dat/TGL_mtl/Longitudinal/).

After that, please execute both ['src/parser_sgl_LS.R'](src/parser_sgl_LS.R) and ['src/parser_mtl_LS.R'](src/parser_mtl_LS.R) files, which will generate a processed data combined two different feature sets (MRI, MRI+META) and two different objective scores (MMSE, ADAS-Cog: TOTAL11).

At the end, please execute both ['src/mergeTP_sgl.py'](src/mergeTP_sgl.py) and ['src/mergeTP_mtl.py'](src/mergeTP_mtl.py) files, which will generate all needed tasks for later model learning.

The number of instances for MRI and MRI+META are shown as two figures below.

1. Number of instances with MRI feature

![](img/MRI_num_instance.png)

2. Number of instances with MRI+META feature

![](img/MRI+META_num_instance.png)



## Experiments and Results
All the experiment is conducted under five-fold cross-validation and are repeated with five different seeds; furthermore, in each fold, experiment results are reweighed by the sample size in each time point.

### Which Algorithms Are Implemented?
For one-dimensional algorithms, we implemented three algorithms as following.
1. Least_Lasso_sgl.m
2. Least_TGL_sgl.m
3. Least_CFGLasso_sgl.m

For two-dimensional algorithms, we implemented five algorithms as follows.
1. Least_Lasso_mtl.m
2. Least_TGL_mtl_g.m
3. Least_TGL_mtl_n.m
4. Least_CFGLasso_mtl_g_fl.m
5. Least_CFGLasso_mtl_n.m

### How To Run Experiments?
For experiments, please refer to ['src/run_Ridge_sgl.m'](src/run_Ridge_sgl.m), ['src/run_Least_sgl.m'](src/run_Least_sgl.m) and ['src/run_Least_mtl.m'](src/run_Least_mtl.m). So far, one still needs to manually comment and uncomment some blocks of codes to examine one specific algorithm. It is noticeable that ['src/run_Ridge_sgl.m'](src/run_Ridge_sgl.m) is used for Ridge algorothm, ['src/run_Least_sgl.m'](src/run_Least_sgl.m) is used for one-dimensional algorithms; while ['src/run_Least_mtl.m'](src/run_Least_mtl.m) is designed for two-dimensional algorithms.



### How To Explain Results?

## Appendix: Experiment Results
1. Correlated Coefficient (CC) Evaluation Metric
![](img/result_coefficient.png)

2. Root Mean Squared Error (RMSE) Evaluation Metric
![](img/result_rmse.png)

## Appendix: Repository Structure
```
.
├── README.md
├── README.md.ipynb
├── dat
│   ├── Origin
│   ├── TGL_mtl
│   │   └── Longitudinal
│   └── TGL_sgl
│       └── Longitudinal
├── img
│   ├── META_feature_info.png
│   ├── MRI+META_num_instance.png
│   ├── MRI_num_instance.png
│   ├── result_coefficient.png
│   └── result_rmse.png
└── src
    ├── Least_CFGLasso_mtl_g_fl.m
    ├── Least_CFGLasso_mtl_n.m
    ├── Least_CFGLasso_sgl.m
    ├── Least_Lasso_mtl.m
    ├── Least_Lasso_sgl.m
    ├── Least_TGL_mtl_g.m
    ├── Least_TGL_mtl_n.m
    ├── Least_TGL_sgl.m
    ├── TDFusedLasso(gcc-old).cpp
    ├── TDFusedLasso.cpp
    ├── TDFusedLasso.mexmaci64
    ├── TDFusedLasso.o
    ├── mergeTP_mtl.py
    ├── mergeTP_sgl.py
    ├── parser_mtl_LS.R
    ├── parser_sgl_LS.R
    ├── run_Least_mtl.m
    ├── run_Least_sgl.m
    ├── run_Ridge_sgl.m
    ├── sll_opts.m
    └── toy-example
        ├── 1d.cpp
        ├── 1d.in
        ├── 2d.cpp
        └── 2d.in

9 directories, 31 files
```

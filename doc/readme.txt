* How to Generate the Dataset?

** I already include required dataset in "report/dat/" folder, including the files for one-dimensional and two-dimensional algorithms. However, it is still possible for you to regenerate all the files with "parser_LS_new" and "parser_mtl_LS_new". The former one generates dataset for one-dimensional algorithms and the latter one is for two-dimensional ones.



* How To Run Experiments?

** For all experiments, they can be categorized into one-dimensional algorithms and two-dimensional algorithms. For the former one, they include 1D-Ridge, 1D-Lasso and 1D-TGL. For the implementation convenience, I just run 1D-Ridge experiments in function "run_Ridge_sgl", and run the rest of them in "run_Least_sgl". On the other hand, 2D-Lasso, 2D-TGL and 2D-TGL+ are implemented in "run_Least_mtl".

** The experiment framework in "run_Ridge_sgl" is simple. There are total of eight experiments consisted by the combination of two feature sets, two target scores, and two evaluation metrics, which produces eight experiment results (2 x 2 x 2).

** The experiment framework in "run_Least_sgl" and "run_Least_mtl" are similar. There are multiple algorithms in each file, which are listed above, and every algorithm will be executed eight times by different combinations of feature sets, target scores, and evaluation metrics. To run one of the algorithms in each file, just counter-comment out the corresponding block of parameter settings and algorithm statement. For example, if you want to execute 1D-TGL algorithm, in "run_Least_sgl" you can just comment out the block in betweeen line 14 to line 19 and line 113, meanwhile, counter-comment out line 24 to line 29 and line 114. Then "run_Least_sgl" will produce the correct experiment result for you.

** Note that all the experiment is conducted under five-fold cross-validation and are repeated with five different seeds; furthermore, in each fold, experiment results are reweighed by the sample size in each time point.

** Table for Referencing algorithms to MatLab functions
*** 1D-Ridge: run_Ridge_sgl
*** 1D-Lasso: run_Least_sgl -> Least_Lasso_sgl
*** 1D-TGL  : run_Least_sgl -> Least_TGL_sgl
*** 1D-cFSGL: run_Least_sgl -> Least_cFSGL_sgl
*** 2D-Lasso: run_Least_mtl -> Least_Lasso_mtl
*** 2D-TGL  : run_Least_mtl -> Least_TGL_mtl_n
*** 2D-TGL+ : run_Least_mtl -> Least_TGL_mtl_g



* How To Explain Experiment Results?

** For each of algorithm experiment, there are eight blocks for different combinations between feature sets, target scores and evaluation metrics. Each of the block contains the task description, parameter description, results for five different seeds, the overall averaged results and their variance.



** How to Examine the Improvement Significance?

** I consider the comparisons between pairs like (1D-Lasso, 2D-Lasso) and (1D-TGL, 2D-TGL) for "Hypothesis I: 2D algorithms defeat 1D algorithms" and pairs like (2D-TGL, 2D-TGL+) for "Hypothesis II: 2D+ algorithms defeats 2D algorithms". All the comparisons are conducted under paired t-test with five samples for both competitors generated from five different seeds (1 to 5).



* How to Understand the Implementation for Each Algorithm?

** Since most of my codes have similar structures, I only add comments on few of them, such as "run_Least_mtl" and "Least_TGL_mtl_g". I believe the rest of them are straightforward to understand.

** For more details, please refer to my previous slides "ITRI_Report" under "report/doc/" and corresponding papers. Or, you can feel free to email your questions to me: "wrangle1005@gmail.com".



* What is the Update Information?

** Organize my project structure 

** Rerun all the experiments with different seeds instead of using merely one seed.

** Show the significance result



* Appendix :: Experiment

** Experiment Results (Correlated Coefficient)
         |     (M, MMSE)     |    (M+E, MMSE)    |     (M, ADAS)     |  (M+E, ADAS)
1D-Ridge | 0.741961±0.000068 | 0.777012±0.000095 | 0.744862±0.000021 | 0.795173±0.000030
1D-Lasso | 0.740918±0.000004 | 0.814807±0.000007 | 0.720986±0.000001 | 0.835312±0.000011
1D-TGL   | 0.778166±0.000029 | 0.831581±0.000004 | 0.769211±0.000022 | 0.850624±0.000014
1D-cFSGL | 0.787235±0.000011 | 0.834423±0.000004 | 0.782664±0.000017 | 0.850632±0.000013
2D-Lasso | 0.755110±0.000019 | 0.826931±0.000018 | 0.730911±0.000006 | 0.834281±0.000002
2D-TGL   | 0.787764±0.000058 | 0.837247±0.000024 | 0.767291±0.000015 | 0.850178±0.000004
2D-TGL+  | 0.817869±0.000002 | 0.865517±0.000001 | 0.801083±0.000002 | 0.873725±0.000002

** Improvement Significance (Correlated Coefficient)
                          | (M, MMSE) |(M+E, MMSE)| (M, ADAS) |(M+E, ADAS)
H1 - (1D-Lasso, 2D-Lasso) | 0.002046  | 0.000235* | 0.002898* |   - - - 
H1 - (1D-TGL  , 2D-TGL  ) | 0.076539  | 0.083429  |   - - -   |   - - - 
H2 - (2D-TGL  , 2D-TGL+ ) | 0.001751* | 0.000213* | 0.000018* | 0.000010* 

** Experiment Results (Root Mean Squared Error)
         |     (M, MMSE)     |    (M+E, MMSE)    |     (M, ADAS)     |  (M+E, ADAS)
1D-Ridge | 3.018537±0.003030 | 2.821958±0.002990 | 6.061067±0.004157 | 5.499035±0.006388
1D-Lasso | 2.952054±0.000137 | 2.532943±0.000044 | 6.142708±0.000546 | 4.855201±0.000929
1D-TGL   | 2.768080±0.001932 | 2.413267±0.000533 | 5.689965±0.004379 | 4.676146±0.000993
1D-cFSGL | 2.739402±0.001347 | 2.412433±0.000366 | 5.638594±0.002056 | 4.595167±0.001067
2D-Lasso | 2.974263±0.000549 | 2.516909±0.000769 | 6.466563±0.000166 | 4.981515±0.000753
2D-TGL   | 2.827151±0.014793 | 2.385045±0.000056 | 5.900847±0.001264 | 4.792991±0.000264
2D-TGL+  | 2.645794±0.001399 | 2.224907±0.000258 | 5.509920±0.001204 | 4.406653±0.000019

** Improvement Significance (Correlated Coefficient)
                          | (M, MMSE) |(M+E, MMSE)| (M, ADAS) |(M+E, ADAS)
H1 - (1D-Lasso, 2D-Lasso) |   - - -   |   - - -   |   - - -   |   - - -  
H1 - (1D-TGL  , 2D-TGL  ) |   - - -   |   - - -   |   - - -   |   - - -  
H2 - (2D-TGL  , 2D-TGL+ ) | 0.022604* | 0.000096* | 0.000016* | <0.000010*

*** Note that "*" means sinificant improvement, while "- - -" means performance decline.
*** Hypothesis I seems work sometimes under Correlated Coefficient evaluation metric; however it fails under Root Mean Squared Error evaluation metric.
*** Hypothesis II seems work everywhere.



* Appendix :: Package Contents
report/
├── MALSAR1.1
│   ├── COPYRIGHT
│   ├── INSTALL.m
│   ├── MALSAR
│   ├── data
│   ├── examples
│   ├── manual
│   └── readme.txt
├── dat
│   ├── Origin
│   ├── TGL
│   └── TGL_mtl
├── doc
│   ├── ITRI_Report.pdf
│   └── readme.txt
└── src
    ├── Least_CFGLasso_sgl.m
    ├── Least_Lasso_mtl.m
    ├── Least_Lasso_sgl.m
    ├── Least_TGL_mtl_g.m
    ├── Least_TGL_mtl_n.m
    ├── Least_TGL_sgl.m
    ├── parser_LS_new.R
    ├── parser_mtl_LS_new.R
    ├── run_Least_mtl.m
    ├── run_Least_sgl.m
    └── run_Ridge_sgl.m

11 directories, 15 files

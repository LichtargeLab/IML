# IML

Integrative Machine Learning (IML)

### Websites

[Lichtarge Lab](http://lichtargelab.org)

[Evolutionary Action](http://eaction/lichtargelab.org)

[Lichtarge Lab Github](https://github.com/lichtargelab)

##### Related publications

[Parvandeh et al., Consensus features nested cross-validation, Bioinformatics, 2020](https://doi.org/10.1093/bioinformatics/btaa046)

[Parvandeh et al., EPIMUTESTR: a nearest neighbor machine learning approach to predict cancer driver genes from the evolutionary action of coding variants, 2022](https://doi.org/10.1093/nar/gkac215)

[Katsonis P, Lichtarge O., A formal perturbation equation between genotype and phenotype determines the evolutionary action of protein coding variations on fitness, Genome Research, 2014 Sep 12. pii: gr.176214.114.](https://pubmed.ncbi.nlm.nih.gov/25217195/)

### To run
1. Download the four data types from [GDSC website](https://cellmodelpassports.sanger.ac.uk/downloads)

2. Annotate mutation data with [Evolutionary Action](http://eaction.lichtargelab.org)

3. Preprocess four data types

      > Rscript 1_Preprocessing_data.R
      
4. Predict IC50 values for pan-cancer cell lines

      > Rscript 2_IML_PANCAN_1Dtype.R
      > Rscript 3_IML_PANCAN_2Dtypes.R
      > Rscript 4_IML_PANCAN_3Dtypes.R
      > Rscript 5_IML_PANCAN_4Dtypes.R
      
5. Predict IC50 values for cancer specific cell lines

      > Rscript 6_IML_CASpecific_1Datatype.R
      > Rscript 7_IML_CASpecific_2Datatypes.R
      > Rscript 8_IML_CASpecific_3Datatypes.R
      > Rscript 9_IML_CASpecific_4Datatypes.R
      
6. Define sensitivity signature using genes that identified from step 5 for RNA-Seq

      > Rscript 10_IML_SensitivitySignature.R
      
7. Visualize the results

      > Rscript 11_Visualization.R
      
8. Validate on CCLE data

      > Rscript 12_Validation_CCLE.R
      
### Dependencies

```
install.packages(c('CORElearn', 'Rcpp', 'dplyr', 'parallel', 'foreach', 'doParallel', 'glmnet', 'randomForest', 'e1071', 'rpart'))

```
Other R packages

```
install.packages(c('ggplot2', 'tidyverse', 'hrbrthemes', 'viridis', 'ggpubr', 'ggbeeswarm', 'forcats', 'cvTools'))

```

Helper function (included)

```
Rcpp::sourceCpp("VCF2DM.cpp")

```

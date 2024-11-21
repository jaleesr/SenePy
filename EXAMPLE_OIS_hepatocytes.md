OIS senePy in R example
================

**NOTE:** It is recommended to use senePy in python for seamless
integration into your workflow. You can save the results as a dataframe
and open it in R.

But if you really want to avoid using python, you can use this workflow
below. However, you will be limited on some functionality.

``` r
library(Seurat)
```

    ## Loading required package: SeuratObject

    ## Loading required package: sp

    ## 
    ## Attaching package: 'SeuratObject'

    ## The following object is masked from 'package:base':
    ## 
    ##     intersect

``` r
library(reticulate)
```

Recommended to install senepy into a conta environment: \$conda create
-n senepy python=3.9 \$conda activate senepy \$pip install senepy

``` r
use_condaenv("senepy", required = TRUE) #name of the conda environment, in this case it is also named senepy

sp <- import("senepy")
```

### Dataset: Chan et al 2024

KRAS oncogene induced cellular senescence in mouse hepatocytes  
10.1038/s41586-024-07797-z

**Samples:** mV: control, D12: 12 days after KRAS induction via
tail-vein injection

Toy data are subset to only include 500 total cells from mV and D12

``` r
seurat_obj <- readRDS("toy_data/TOY_DATA_Chan_et_al_2024.rds")
seurat_obj
```

    ## An object of class Seurat 
    ## 15766 features across 500 samples within 1 assay 
    ## Active assay: RNA (15766 features, 0 variable features)
    ##  1 layer present: counts

``` r
hubs <- sp$load_hubs(species = 'Mouse')
```

### Deciding which signature to use

These are mouse hepatocytes, so we want to focus on liver/hepatocytes

``` r
unique(hubs$metadata$tissue) #to get all unique tissues
```

    ##  [1] "Bladder"         "Brain"           "Diaphragm"       "Fat"            
    ##  [5] "HSC"             "Heart_and_Aorta" "Kidney"          "Large_Intestine"
    ##  [9] "Limb_Muscle"     "Liver"           "Lung"            "Lymphoid"       
    ## [13] "Myeloid"         "Pancreas"        "Thymus"          "Tongue"         
    ## [17] "Trachea"

``` r
hubs$metadata[hubs$metadata$tissue == 'Liver',]
```

    ##    tissue                                 cell hub_num size n_sen          hyp
    ## 41  Liver                         Kupffer cell       0  638    14 2.268194e-03
    ## 42  Liver                         Kupffer cell       1 3305    17 9.983103e-01
    ## 43  Liver endothelial cell of hepatic sinusoid       0  205    12 4.280991e-07
    ## 44  Liver                           hepatocyte       0 1600    14 6.086587e-01
    ## 45  Liver                           hepatocyte       1  143     7 3.511115e-04

**Making the unviversal signature**

We pass hubs.metadata, which contains every senePy signature

We specify **calculate_thresh = True** to find the genes that are
statistically overrepresented and set the FDR p-value threshold to 0.01

``` r
merge_results <- hubs$merge_hubs(hubs$metadata,
                                 new_name = 'Universal',
                                 calculate_thresh = TRUE,
                                 p_thres = 0.01)
```

    ## A gene will occur 5 times at 0.44% chance
    ## Threfore 5 is the calculated_threshold

to get a list of genes based on a signature:

``` r
universal_sig <- sapply(hubs$hubs[["Universal"]], function(x) x[[1]])
```

if we wanted a cell-specific signature we need to pass the whole python
tuple as a string

``` r
hep_1 <- hubs$hubs[["('Liver', 'hepatocyte', 1)"]]
hep_1 <-sapply(hep_1, function(x) x[[1]]) #if you can use the weights then you can grab them with x[[2]]
```

We have extracted signatures of interest but we cannot use senePy
directly against a Seurat object. So we will have to use the Seurat
module score.

Note, that we also cannot use the weights, which are an important
compenent of the senepY signature when scoring cells.

``` r
seurat_obj <- NormalizeData(seurat_obj)
```

    ## Normalizing layer: counts

``` r
seurat_obj <- AddModuleScore(
  object = seurat_obj,
  features = list(universal_sig),
  name = "universal_score"
)
```

    ## Warning: The following features are not present in the object: Hoxc4, Ablim3,
    ## Fgf9, Slc38a5, 1700019D03Rik, 4931408D14Rik, Ccl11, Lrrc26, Ctsg, H2-Q8,
    ## Clec7a, Gm7609, Igj, Gm11428, Atp6v0c-ps2, 1810033B17Rik, Lilrb4,
    ## 2010001M09Rik, Chi3l1, Stmn2, Cd209f, Tmem40, AF251705, 6030408B16Rik,
    ## 1100001G20Rik, Fcnb, Gp49a, Lilrb3, Itgb2l, A430084P05Rik, Olfr613, Folr4,
    ## Ms4a3, Sfpi1, Ankrd22, Hmha1, Gm4902, 6330512M04Rik, C7, Htra4, Wisp2, Mnda,
    ## Sostdc1, Gm12250, Faim3, Amica1, Retnla, H2-T10, Gm12504, Hbb-b1, Hbb-b2,
    ## BC064078, Gm5069, 4930572J05Rik, 4930506M07Rik, Kdm5d, Mpl, Uchl1, C1qtnf9,
    ## Sncg, Rbp7, Timp4, Gypa, Krtdap, Gm13315, Fam151a, Resp18, Them5, Prss16,
    ## Ccdc109b, Tmem27, Slc4a1, Calml3, Reg3g, Rorb, Camp, Sftpa1, Ms4a4d, Chi3l3,
    ## Sftpc, Krt84, Mpo, Krt36, Elane, Car4, Gpr116, Tcf15, BC020535, Aqp2, BC021785,
    ## Mx2, 1700112E06Rik, 4930413G21Rik, Ngp, Ltf, Mmp3, Fam71a, Cbr2, Slc13a1,
    ## Cyp2d12, Rps19-ps3, Tmigd1, Slc22a12, Nat8, Zpbp, D730005E14Rik, Upk1a,
    ## S100a14, Oas1b, Dem1, Niacr1, Aqp3, 1190002F15Rik, Beta-s, Aldh3b2, Slfn10-ps,
    ## Fam26f, Srpk3, Rnls, 0610010O12Rik, Prss34, Ptx3, Nppa, Actc1, not searching
    ## for symbol synonyms

``` r
#note that for some reason AddModuleScore adds a number to the end of the score name, i.e., the resulting score name is universal_score1
```

``` r
#adding a condition column to the seurat object metadata based on sample name
seurat_obj@meta.data$Condition <- sub("_.*", "", seurat_obj@meta.data$Sample)
```

``` r
wilcox.test(universal_score1 ~ Condition, data = seurat_obj@meta.data)
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  universal_score1 by Condition
    ## W = 31299, p-value = 0.02542
    ## alternative hypothesis: true location shift is not equal to 0

We still get a significant shift based on this gene set, but it is not
as pronounced when the genes are not weighted by their importance. Also,
their is no gene name-harmonization between senePy and the object.
Agian, we recommend using senePy directly in python.

Example of gene set enrichment

``` r
Idents(seurat_obj) <- "Condition"
seurat_obj$Condition <- as.factor(seurat_obj$Condition)

#not actually a great way to do DE in single-cell analysis... just for example purporses
de_results <- FindMarkers(
  object = seurat_obj,
  ident.1 = "D12",
  ident.2 = "mV",
  test.use = "wilcox"
)
```

    ## For a (much!) faster implementation of the Wilcoxon Rank Sum Test,
    ## (default method for FindMarkers) please install the presto package
    ## --------------------------------------------
    ## install.packages('devtools')
    ## devtools::install_github('immunogenomics/presto')
    ## --------------------------------------------
    ## After installation of presto, Seurat will automatically use the more 
    ## efficient implementation (no further action necessary).
    ## This message will be shown once per session

``` r
de_results <- de_results[de_results$p_val_adj < 0.05 & de_results$avg_log2FC > 0,]

head(de_results)
```

    ##               p_val avg_log2FC pct.1 pct.2    p_val_adj
    ## Lcp1   1.930109e-19   3.244268 0.416 0.018 3.043009e-15
    ## Scand1 3.197282e-19   1.550691 0.596 0.125 5.040835e-15
    ## Slfn2  5.359582e-19   3.038953 0.416 0.024 8.449917e-15
    ## Rbm3   1.170846e-18   1.396087 0.560 0.107 1.845956e-14
    ## Ifi30  2.430208e-18   3.471844 0.413 0.030 3.831466e-14
    ## Psmb9  2.850508e-18   1.906785 0.461 0.048 4.494111e-14

Simple hypergeometic gene-set enrichment. You can use the the signature
directly in any enrichment method you want.

``` r
background_genes <- rownames(seurat_obj[["RNA"]])
N <- length(background_genes)

gene_set_in_background <- intersect(universal_sig, background_genes)
K <- length(gene_set_in_background)

de_genes <- rownames(de_results)
n <- length(de_genes)

k <- length(intersect(de_genes, gene_set_in_background))
```

``` r
phyper(q = k - 1, m = K, n = N - K, k = n, lower.tail = FALSE)
```

    ## [1] 0.006068899

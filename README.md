# [Identification of multivariate phenotypes most influenced by mutation: _Drosophila serrata_ wings as a case study](https://doi.org/10.1093/evolut/qpaf160) <a href="http://rsbl.royalsocietypublishing.org/content/12/5/20151003"><img src="http://tguillerme.github.io/images/OA.png" height="20" widht="20"/></a>
[Cara Conradsen](https://mcguiganlab.org/cara-conradson/) and [Katrina McGuigan](https://mcguiganlab.org/katrina-mcguigan/)

Published in Evolution: 01 August 2025

This repository contains all the code and data used in the manuscript.

###### We are still in the process of finalising and refining the manuscript.
###### The history of the experiment can be found in [Causes of variability in estimates of mutational variance from mutation accumulation experiments](https://doi.org/10.1093/genetics/iyac060) in Genetics. <a href="http://rsbl.royalsocietypublishing.org/content/12/5/20151003"><img src="http://tguillerme.github.io/images/OA.png" height="20" widht="20"/></a>

## Analysis
### `MVN_M_6_UN.sas`
This SAS script analyses the Drosophila wing trait data from our mutation reduction experiment where the 12 populations were analysed separately. The outcome is the multivariate (MVN) unstructured (UN) mutation variance - covariance matrix (M) for six wing traits. Here, populations are labelled first by treatment and then sampled generation (e.g., treatment 1_1 small treatment, generation 1). The steps include:

1. Sorting the Data
The dataset `mrwings` is sorted by `Trait` and `Gen`.
2. Standardising Scores
The `Score` variable is standardised (mean = 0, std = 1) within each `Trait` group and saved as `stdmrwings`.
3. Linear Mixed Model using the mixed procedure (`PROC MIXED`)
A mixed-effects model is fit using REML, applied to treatments seperately, where multivaraite outliers are excluded (e.g., `Treat=1`, `Gen=1`, excluding `multiout=1`):
- _Fixed effect: `Trait`
- _Random effects: `Trait` within `Line`, `Trait` within `Vial`, nested in `Treat` and `Line`
- _Repeated measures: `Animal` nested in `Treat`, `Line`, and `Vial`
4. Output Tables
The following outputs are saved for downstream analysis (within the folder `un_sas_output`):
- _`Covparms1_1`: Covariance parameters
- _`FitStat1_1`: Fit statistics
- _`Converge1_1`: Convergence diagnostics
- _`Asycov1_1`: Asymptotic covariance matrix

### `randomise_wing_data.R` 

There are four parts to `randomise_wing_data.R`. In the first part, _Part 1. Prepare the data set to be shuffled_, we create a unique ID for a line. This section also separates the experimental component of Generation + Treatment + Line (`population_info`) and the individual fly trait information (with a newly assigned unique line ID, `line_trait_info`) into two data frames. 

In the second section (_Part 2. Generate a data frame of sets of randomised unique line IDs_), here are the steps where we randomly sample the list of unique line IDs, and then generate a data frame of the sampled line sets (`all_sets`, _N_ = 1000), each column is a single set of sampled lines and each row, within a column, is a randomly sampled line. This is where you can set the number of data sets needed and whether you want to sample with or without replacement. 

In the third part (_Part 3. Write a function to merge randomised line data to trait info and save to .csv_), we create the function `create_randomise_wing_data` that uses the data frames generated in parts 1 and 2. It takes a single column from `all_sets` (one randomised sample of lines),  and does several things:
1. adds the randomised unique line IDs by pasting the `all_sets` column onto the data frame of the experimental component, `population_info`
- _The unique lines are now randomised across the 12 experimental populations of treatment X generation and 42 experimental lines._
2. then populates with the flies and their wing information by appending the assigned unique line ID in `line_trait_info` to the randomised unique line ID
- _This keeps replicate vials and individual flies intact within line, and at this step, this is our complete randomised dataset_
3. we retain the column with the experimental design lines, numbered from 1 to 42, for SAS and remove the column with the unique line IDs, numbered from 1 to 493 (_see Note below_)
4. calculates the Mahalanobis Distance and determines multivariate outliers for this new data set to use in SAS
5. adds a column of the treatments as numbers for SAS, using the original notation: small is 1 and big is 3 (SB crosses were 2)
6. converts the new data set into long format for SAS
7. re-orders the data columns
8. sorts the data by Gen, Treatment, Line, Vial, Individual
9. saves the dataset as a .csv in the 'generated_datasets' subfolder

Finally, in the last section (_Part 4. Implement function on sets of randomised lines_) we implement the function created above on the data frame `all_sets`. In this section, to speed up the time to generate all 1000 data sets, the function is run in parallel using the `foreach` function (%dopar%).  

> [!Note]
> Our dataset had fewer 'unique' lines than the full experimental design (i.e., 6 gens X 2 treat X 42 lines = 504, Observed data = 493). Our code is written to keep the same number of 'unique' lines (_N_ = 493), where shuffled genetic data is extended onto this.
> If sample with replacement equals:
> - False; same number of individuals (_N_ = 5135)
> - True; the number of total individuals will vary

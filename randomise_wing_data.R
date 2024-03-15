#### RANDOMISE WING DATA SETS

### Function to generate N data sets to sample from to generate a Null distribution.
### We would generate ~1,000 data sets where we randomly sampled (with or without 
### replacement) a line (both reps and all individuals per) and assigned it to one 
### of 12 groups (6 generations, 2 treatments), and then estimated the M for each of 
### the 12 groups.

# [1] "Animal"    "Gen"       "Treatment" "Treat"     "Line"      "Vial"      "Ind"       "MD"       
# [9] "multiout"  "Index"     "Trait"     "Score" 

# Packages
library(data.table); library(foreach); library(parallel); library(doParallel)

# create an output directory
outdir <- "generated_datasets"

if (file.exists(outdir)==FALSE){
  dir.create(outdir)
}

# Import wing data
wing_dat <- fread("mr_wings_6traits.csv")

# Create a number key for traits
trait_number <- unique(wing_dat[,c("Index", "Trait")])

### Part 1. Prepare the data set to be shuffled ####
# Outcome will be two data tables
# population_info: Gen, Treatment, Line
# line_trait_info: Line +vial +ind and Score

#remove unnecessary columns
rm_col <- c("stdScore","outSD","multiout", "MD", "Trait")
set(wing_dat, ,rm_col, NULL)

# convert to wide
wing_data_wide <- dcast(wing_dat, Animal+Gen+Treatment+Treat+Line+Vial+Ind ~ Index,  
                        value.var ="Score")

# Population Information
# Get the information for Gen, Treatment and Line
population_info <- unique(wing_data_wide[, c("Gen", "Treatment", "Line")])
nrow(population_info) # 493 unique 'lines'

# Because line 42 in gen 1 S is different from line 42 in gen 2 B, 
# need to create a unique line identifier which we will use to randomise lines
population_info$LineID = 1:nrow(population_info)

# Add the unique line identifier back to the observed data, 
# omitting missing lines
wing_data_wide <- merge(wing_data_wide, population_info, 
                        by=c("Gen","Treatment", "Line"), all.x = TRUE)

# check number of unique lines
length(unique(wing_data_wide$LineID)) # should be 493

  
# Line and Trait Information
line_trait_info <- wing_data_wide[,!c("Gen","Treatment","Treat", "Line")]# 

# Remove LineID from population_info (this will be ran)
population_info <- population_info[, .SD, .SDcols = !"LineID"]

### Part 2. Write a function to randomised data and save to .csv ####

# pop_dat=population_info; biol_dat=line_trait_info; line_smpl=all_sets[,1]


create_randomise_wing_data <- function(line_smpl,
                                pop_dat=population_info, biol_dat=line_trait_info){
  # set default pop_dat as population_info and  biol_dat as line_trait_info
  rep = names(line_smpl)
  pop_dat$LineID = line_smpl
  randomised_line_data <- biol_dat[pop_dat, on = "LineID"]
  randomised_line_data$LineID <- NULL # remove LineID
  
  # Add Mahalanobis Distance and outliers
  x <- randomised_line_data[, .SD, .SDcols = c("CS","ILD1.2","ILD1.5","ILD2.5","ILD2.8","ILD3.7")]
  randomised_line_data$MD <- mahalanobis(x, colMeans(x), cov(x))
  chi_crit_val <- qchisq(0.999, df=6)
  randomised_line_data$multiout <- randomised_line_data[,fcase(MD >= chi_crit_val, 1,
                                                               MD < chi_crit_val, 0, 
                                                               default = NA)]
  
  # Add treatment as an integer "S" = 1 and "B" = 3
  randomised_line_data$Treat <- randomised_line_data[,fcase(Treatment=="S", 1,
                                                            Treatment=="B", 3, 
                                                               default = NA)]
  
  # Convert to long
  randomised_line_data_lng <- melt(randomised_line_data, 
                                   id.vars=c("Gen","Treatment","Treat","Line", 
                                             "Animal", "Vial", "Ind", "MD", "multiout"),
                                   variable.name = "Index", value.name = "Score")
  
  # Add integer code for traits
  randomised_line_data_lng <- trait_number[randomised_line_data_lng, on = "Index"]
  

  # re-order columns
  setcolorder(randomised_line_data_lng,
           c("Animal", "Gen","Treatment", "Treat", "Line",  "Vial","Ind", "MD", 
             "multiout", "Index", "Trait", "Score"))
  
  # sort
  setorderv(randomised_line_data_lng,
            c("Gen", "Treat", "Line","Vial", "Ind"))
  
  # Output csv
  fwrite(randomised_line_data_lng, 
         file=paste0(outdir, "/mr_wings_6traits_rando_", rep, ".csv"),
         quote=FALSE)

}


### Implement function on sets of randomised lines ####

# Create data frame of randomised sets of 'lines'

# Set number of data sets to generate
n_set = 20
# sample with replacement?
smpl_with_replace = TRUE

# Create a data.frame of random line samples,
# where each column is a unique draw of line ids 1:493
line_id_list <- unique(line_trait_info$LineID)
num_lines = length(line_id_list)

all_sets <-data.table(NULL)
for (i in 1:n_set) {
  sub_set <- sample(line_id_list, size = num_lines, replace = smpl_with_replace)
  all_sets <- cbind(all_sets,sub_set)
}
colnames(all_sets) <- as.character(seq(1,n_set))

# Detect and set cores for parallel processing
n_cores <- detectCores()-1
cluster <- makeCluster(n_cores)
doParallel::registerDoParallel(cluster)

# Implement function in parallel
foreach(x = 1:ncol(all_sets), .packages = c("data.table")) %dopar%
  create_randomise_wing_data(all_sets[,..x])


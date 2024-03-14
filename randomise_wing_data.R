#### RANDOMISE WING DATA SETS

### Function to generate N data sets to sample from to generate a Null distribution.
### We would generate ~1,000 data sets where we randomly sampled (with or without 
### replacement) a line (both reps and all individuals per) and assigned it to one 
### of 12 groups (6 generations, 2 treatments), and then estimated the M for each of 
### the 12 groups.

# [1] "Animal"    "Gen"       "Treatment" "Treat"     "Line"      "Vial"      "Ind"       "MD"       
# [9] "multiout"  "Index"     "Trait"     "Score"     "stdScore"  "outSD" 

# Packages
library(data.table); library(purrr); library(parallel); library(dplyr)

# create an output directory
dir.create("generated_datasets")

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
# need to create a unique line identifier
population_info$LineID = 1:nrow(population_info)

# Add the unique line identifier back to the observed data, 
# omitting missing lines
wing_data_wide <- merge(wing_data_wide, population_info, 
                        by=c("Gen","Treatment", "Line"), all.x = TRUE)

# check number of unique lines
length(unique(wing_data_wide$LineID)) # should be 493

  
# Line and Trait Information
line_trait_info <- wing_data_wide[,!c("Gen","Treatment","Treat", "Line")]

# Remove LineID from population_info (this will be ran)
population_info <- population_info[, .SD, .SDcols = !"LineID"]

### Part 2. Randomise datasets ####

# start_time <- Sys.time()
# Set number of data sets to generate
n_set = 5

# Create a data.frame of random line samples,
# where each column is a unique draw of line ids 1:493
line_id_list <- unique(line_trait_info$LineID)
num_lines = length(line_id_list)

all_sets <-data.table(NULL)
for (i in 1:n_set) {
  sub_set <- sample(line_id_list, size = num_lines, replace = FALSE)
  all_sets <- cbind(all_sets,sub_set)
}
colnames(all_sets) <- as.character(seq(1,n_set))
# end_time <- Sys.time()
# end_time - start_time
# 1000 number of sets takes ~7.073825 secs

# Determine number of uniquely sampled 'lines'
#length(unique(unlist(c(all_sets[,1]))))

create_randomise_wing_data <- function(line_smpl,
                                pop_dat=population_info, biol_dat=line_trait_info){
  # set default pop_dat as population_info and  biol_dat as line_trait_info
  rep = names(line_smpl)
  pop_dat$LineID = line_smpl
  merge(pop_dat, biol_dat, by="LineID", all.y = TRUE)
}



# Check for multivariate outliers
# MD multiout
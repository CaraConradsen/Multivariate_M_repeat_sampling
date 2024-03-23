#### RANDOMISE WING DATA SETS

### Function to generate N data sets to sample from to generate a Null distribution.
### We would generate ~1,000 data sets where we randomly sampled (with or without 
### replacement) a line (both reps and all individuals per) and assigned it to one 
### of 12 groups (6 generations, 2 treatments), and then estimated the M for each of 
### the 12 groups.

# Packages
library(data.table); library(foreach); library(parallel); library(doParallel)

# create an output directory
outdir <- "generated_dataset"

if (file.exists(outdir)==FALSE){
  dir.create(outdir)
}

### Part 1. Prepare the data set to be shuffled ####
# Outcome will be two data tables
# population_info: Gen, Treatment, (observed) Line
# line_trait_info: animal +vial +ind and Trait score in wide format and a unique LineID

# Import wing data
wing_dat <- fread("mr_wings_6traits.csv")

# Create a number key for traits
trait_number <- unique(wing_dat[,c("Index", "Trait")])

#remove unnecessary columns
rm_col <- c("stdScore","outSD","multiout", "MD", "Trait")
wing_dat <- wing_dat[, .SD, .SDcols = !rm_col]

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
# omitting lines from the observed data that had no observations
wing_data_wide <- merge(wing_data_wide, population_info, 
                        by=c("Gen","Treatment", "Line"), all.x = TRUE)

# check number of unique lines
length(unique(wing_data_wide$LineID)) # should be 493

  
# Line and Trait Information
line_trait_info <- wing_data_wide[,!c("Gen","Treatment","Treat", "Line")]# 

# Remove LineID from population_info (this will be ran)
population_info <- population_info[, .SD, .SDcols = !"LineID"]

### Part 2. Generate a data frame of sets of randomised unique line IDs ####
# Create data frame of randomised sets of 'lines' using the unique line ID
# all_sets: each column is a sample of lines [set] and each row is a line

# Set number of data sets to generate
n_set = 1000
# sample with replacement? TRUE/FALSE
smpl_with_replace = TRUE
# set seed for reproducibility
set.seed(42)

# Create a list of the unique line IDs
line_id_list <- unique(line_trait_info$LineID)
# set the number of lines to be sampled
num_lines = length(line_id_list)# here 493

# Create a data.frame of random line samples,
# where each column is a unique draw of line ids 1:493

all_sets <-data.table(NULL)
for (i in 1:n_set) {
  sub_set <- sample(line_id_list, size = num_lines, replace = smpl_with_replace)
  all_sets <- cbind(all_sets,sub_set)
}
# Rename columns to the set number
colnames(all_sets) <- as.character(seq(1,n_set))

# Inspect first 5 columns of all_sets
all_sets[,1:5]

### Part 3. Write a function to merge randomised line data to trait info ####

create_randomise_wing_data <- function(line_smpl,
                                pop_dat=population_info, biol_dat=line_trait_info){
  # set default pop_dat as population_info and  biol_dat as line_trait_info
  pop_dat$RandRep = as.integer(names(line_smpl))# stores the column numbered name to number data set later
  pop_dat$LineID = line_smpl# adds the randomised unique line IDs to population_info
  
  # Populate with the flies and their wing information by merging on "LineID"
  randomised_line_data <- biol_dat[pop_dat, on = "LineID"]
  randomised_line_data$LineID <- NULL # remove the unique LineID, no longer needed
  
  # Add Mahalanobis Distance and outliers
  x <- randomised_line_data[, .SD, .SDcols = c("CS","ILD1.2","ILD1.5","ILD2.5","ILD2.8","ILD3.7")]
  randomised_line_data$MD <- mahalanobis(x, colMeans(x), cov(x))
  chi_crit_val <- qchisq(0.999, df=6)
  randomised_line_data$multiout <- randomised_line_data[,fcase(MD >= chi_crit_val, 1,
                                                               MD < chi_crit_val, 0, 
                                                               default = NA)]
  
  # Add treatment as a number for SAS: "S" = 1 and "B" = 3
  randomised_line_data$Treat <- randomised_line_data[,fcase(Treatment=="S", 1,
                                                            Treatment=="B", 3, 
                                                               default = NA)]
  
  # Convert to long
  randomised_line_data_lng <- melt(randomised_line_data, 
                                   id.vars=c("RandRep","Gen","Treatment","Treat","Line", 
                                             "Animal", "Vial", "Ind", "MD", "multiout"),
                                   variable.name = "Index", value.name = "Score")
  
  # Add integer code for traits
  randomised_line_data_lng <- trait_number[randomised_line_data_lng, on = "Index"]
  

  # re-order columns
  setcolorder(randomised_line_data_lng,
           c("RandRep","Animal", "Gen","Treatment", "Treat", "Line",  "Vial","Ind", "MD", 
             "multiout", "Index", "Trait", "Score"))
}

### Part 4. Implement function on sets of randomised lines ####

# Detect and set cores for parallel processing
n_cores <- detectCores()-1
cluster <- makeCluster(n_cores)
doParallel::registerDoParallel(cluster)


# Implement function in parallel
randomised_wing_dataset<- foreach(x = 1:ncol(all_sets), .packages = c("data.table"), 
                                  .combine = "rbind") %dopar%
  create_randomise_wing_data(all_sets[,..x])

# sort
setorderv(randomised_wing_dataset,
          c("RandRep","Gen", "Treat", "Line","Vial", "Ind"))

# Output csv
fwrite(randomised_wing_dataset,
       file=paste0(outdir, "/mr_wings_6traits_rando_", n_set, "_datasets.csv"),
       quote=FALSE)


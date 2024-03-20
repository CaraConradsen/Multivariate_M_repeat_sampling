# About 

With Part 1. Prepare the data set to be shuffled, which creates the unique ID for a line. This section also separates the experimental component of Generation + Treatment + Line ('population_info') and the individual fly trait information (with a newly assigned unique line ID, 'line_trait_info') into two data frames. 

In the second section (Part 2. Generate a data frame of sets of randomised unique line IDs, line 65), here are the steps where we randomly sample the list of unique line IDs, and then generate a data frame of the sampled line sets ('all_sets', N=1000), each column is a single set of sampled lines and each row, within a column, is a randomly sampled line. This is where you can set the number of data sets needed and whether you want to sample with or without replacement (you said with, but seemed undecided when we spoke last). 

In the third part (Part 3. Write a function to merge randomised line data to trait info and save to .csv, line 93), we create the function 'create_randomise_wing_data' that uses the data frames generated in parts 1 and 2. It takes a single column from 'all_sets' (one randomised sample of lines),  and does several things:
adds the randomised unique line IDs by pasting the 'all_sets' column onto the data frame of the experimental component, 'population_info'
The unique lines are now randomised across the 12 experimental populations of treatment X generation and 42 experimental lines.
then populates with the flies and their wing information by appending the assigned unique line ID in 'line_trait_info' to the randomised unique line ID
This keeps replicate vials and individual flies intact within line, and at this step, this is our complete randomised dataset
we retain the column with the experimental lines, numbered from 1 to 42, for SAS and remove the column with the unique line IDs, numbered from 1 to 493
calculates the Mahalanobis Distance and determines multivariate outliers for this new data set to use in SAS
adds a column of the treatments as numbers for SAS, using the original notation: small is 1 and big is 3 (SB crosses where 2)
converts the new data set into long format for SAS
re-orders the data columns
sorts the data by Gen, Treatment, Line, Vial, Individual
saves the dataset as a .csv in the 'generated_datasets' subfolder

Finally, in the last section ('Part 4. Implement function on sets of randomised lines', line 144) we implement the function created above on the data frame 'all_sets'. In this section, to speed up the time to generate all 1000 data sets, the function is run in parallel using the foreach function (%dopar%).  

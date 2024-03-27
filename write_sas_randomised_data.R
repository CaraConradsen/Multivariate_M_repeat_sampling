# The following code both separates the 1000 randomised data sets into eight chunks, and writes the
# corresponding proc mixed models to analyse the in the desktop version of SAS. 

library(data.table)
outdir <- "generated_dataset"
fread(paste0(outdir,"/mr_wings_6traits_rando_1000_datasets.csv"))

#### Write 8 lots of code to run in SAS ####
eight_sets <- data.frame(minRep=seq(from=0, to=875, length=8)+1,
                         maxRep=seq(from=125, to=1000, length=8))

for (i in 1:8) {
  reps_nums <- seq(eight_sets[i,1],eight_sets[i,2])
  
  # Output csv outdir, 
  fwrite(randomised_wing_dataset[RandRep %chin% reps_nums],
         file=paste0("R:\\Multivariate Dataset Sampler 2024\\mr_wings_6traits_rando_", eight_sets[i,1],"_",eight_sets[i,2], "_datasets.csv"),
         quote=FALSE)
  sink(paste0("R:\\Multivariate Dataset Sampler 2024\\RND_data_M_",eight_sets[i,1],"_",eight_sets[i,2],".txt"))
  cat(
    'proc import datafile="R:\\Multivariate Dataset Sampler 2024\\mr_wings_6traits_rando_',eight_sets[i,1],"_",eight_sets[i,2],'_datasets.csv"
out=RNDmrwings dbms=csv replace; 
getnames=yes; 
run;
\n
\n
proc sort data = RNDmrwings;
by RandRep Trait;
run;
quit;
\n
\n
proc standard data=RNDmrwings out=stdRND std=1 mean=0;
by RandRep Trait;
var Score;
run;
quit;
\n
\n
proc sort data = stdRND;
by RandRep Treat Gen;
run;
quit;
\n
\n
ods output convergencestatus=converge covparms=M;
proc mixed data=stdRND method=reml lognote scoring=5;
where multiout not in (1);
by RandRep Treat Gen;
class Treat Line Animal Vial Trait;
model Score = Trait;
random Trait / subject= Line type=un;
random Trait / subject= Vial(Treat Line) type=un;
repeated / subject = Animal(Treat Line Vial) type=un;
run;
quit;
\n
\n
PROC EXPORT DATA= WORK.M
OUTFILE= "R:\\Multivariate Dataset Sampler 2024\\Randomised data\\covpars_',eight_sets[i,1],"_",eight_sets[i,2],'.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;											
run;											
quit;
PROC EXPORT DATA= WORK.converge
OUTFILE= "R:\\Multivariate Dataset Sampler 2024\\Randomised data\\converge_',eight_sets[i,1],"_",eight_sets[i,2],'.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;											
run;											
quit;', sep = "")
  sink()
  
}

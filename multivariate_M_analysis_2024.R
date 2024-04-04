# Multivariate Analyses of Mutation

#packages
library(dplyr, warn.conflicts = FALSE); options(dplyr.summarise.inform = FALSE); library(data.table)
library(stringi); library(stringr); library(magrittr); library(evolqg); library(Matrix)
library(foreach); library(matrixStats); library(MASS); library(parallel); library(matrixcalc)
library(gdata); library(bigstatsr); library(abind); library(psych); library(stats)

# library(devEMF)
# library(tidyr); library(plyr); library(Rcpp)
# library(RColorBrewer); library(plotrix) 
# library(gridExtra);library(grid);library(xtable);library(errors);;library(lmtest)
# ;library(viridis)#Rank-based fitting of linear models
# library(corpcor); library(svglite)
# fonts <- list(`Times New Roman` = "DejaVu Serif")

# figure outputs
library(sysfonts); library(showtextdb); library(showtext)
## Load Times New Roman fonts on Windows
font_add("times", regular = "times.ttf", bold = "timesbd.ttf", italic = "timesi.ttf", bolditalic ="timesbi.ttf")

# create an output directory
outdir_tab <- "output_tables"
outdir_fig <- "output_figures"

if (file.exists(outdir_tab)==FALSE){
  dir.create(outdir_tab)
}
if (file.exists(outdir_fig)==FALSE){
  dir.create(outdir_fig)
}

# Section 1. Define overall parameters, import data and  check model convergence and output ------------------------------------------------------
# This section will import and inspect the output of Equation 1 implemented in a 
# restricted-maximum likelihood (REML) framework using in PROC MIXED in SAS v9.4.
# At the end of this section we will have the following arrays:
# 1. p observed M; M_array, dim = c(n, n, p)
# 2. p asymptotic variance-covariance matrices; V_array, dim=c((n*(n+1)/2), (n*(n+1)/2), p)
# 3. p theta vectors of observed covariance parameter estimates, theta_array, dim=c((n*(n+1)/2), 1, p)
# 4. 10,000 M matrices generated via REML-MVN, N ~ (theta, V); AsycovM_array, dim=c(n, n, MVNsample, p)
# 5. The M matrices from the 1000 randomised datasets; null_M_array, dim = c(n,n,nullnumber,p), where

# Overall parameters 
n <- 6  # number of traits 
p <- 12  # number of matrices to compare (2 treatments by 6 generations)
traitnumber <- c(1,2,3,5,6,7) # unique trait number
nullnumber <- 750 # number of randomised datasets
MVNsample <- 10000  # number of REML-MVN samples
epsilon = 0.05 # used for plotting caps CI intervals 


# Specify subfolders
un_cov_dir <- list.files(pattern = "un_sas_output", full.names = TRUE)
rando_un_cov_dir <- list.files(pattern = "generated_dataset", full.names = TRUE)

# Check model convergence for observed M and the 1000 randomized M 
# Observed data:
converge_obs <- list.files(path = un_cov_dir, pattern = "Converge")

un_obs_converge <- foreach(i = 1:length(converge_obs), .combine = rbind) %do% 
  fread(paste(un_cov_dir,converge_obs[i], sep = "/"))[, pop_num := converge_obs[i]]

un_obs_converge # It all looks good! :)

# randomised data:
converge_null <- list.files(path = rando_un_cov_dir, pattern ="converge")

un_null_converge <- foreach(i = 1:length(converge_null), .combine = rbind) %do% 
  fread(paste(rando_un_cov_dir,converge_null[i], sep="/"))

un_null_converge[Status == 1]# only five instances where models did not converge

unconverged_null <- un_null_converge[Status == 1]# will omit the unconverged models

# Import the asymptotic variance-covariance matrix, V, from REML
# Because we are looking at a 6x6 unstructured matrix, we get 21 unique covariance parameter estimates, 
# thus the V matrix becomes 21 x 21
Vfilenames <- list.files(path = un_cov_dir, pattern = "Asycov")

V_list <- lapply(Vfilenames, function (x) 
  as.matrix(fread(paste(un_cov_dir, x, sep = "/"))[c(1:21), c(3:23)]))

names(V_list) <- gsub(".csv","", gsub("Asycov","", Vfilenames)) 

# test whether the asymptotic variance-covariance matrices are positive definite 
lapply(V_list, is.positive.definite)
# Again, all looks good, all are positive-definite! :)


# Get observed M matrices from REML
Mfilenames <- list.files(path = un_cov_dir, pattern = "Covparms")

# create a function to import M and convert the column list into a 6 x 6 matrix
import_M <- function(mlist_name){
  Mdf <- fread(paste(un_cov_dir,mlist_name, sep = "/"))[Subject == "Line",.(Estimate)] 
  M <- matrix(NA, nrow = n , ncol = n)
  upperTriangle(M, diag = TRUE) <- Mdf$Estimate
  lowerTriangle(M) = upperTriangle(M, byrow = TRUE)
  return(M)
}
M_list <- lapply(Mfilenames, import_M)

names(M_list) <- gsub(".csv","", gsub("Covparms","", Mfilenames)) 


#create the M and V arrays
M_array <- array(NA, dim = c(n, n, p))
for (i in 1:p){M_array[,,i] <- M_list[[i]]}
dimnames(M_array) <- list(traitnumber, traitnumber, names(M_list))

V_array <- array(NA, dim=c((n*(n+1)/2), (n*(n+1)/2), p))
for (i in 1:p){V_array[,,i] <- V_list[[i]]}
dimnames(V_array)[[3]] <- names(V_list)


# vectors of covparms from REML that will generate REML-MVN estimates of M 
theta_list <- sapply(Mfilenames, function (x)
  fread(paste(un_cov_dir, x, sep = "/"))[Subject == "Line",.(Estimate)])
names(theta_list) <- gsub(".csv", "", gsub("Covparms", "", Mfilenames))  

theta_array <- array(NA, dim=c((n*(n+1)/2), 1, p))
for (i in c(1:p)){theta_array[,,i] <- theta_list[[i]]}
dimnames(theta_array)[[3]] <- names(theta_list)

# Create array for 10,000 REML-MVN data
AsycovM_array<-array(NA, dim=c(n, n, MVNsample, p)) 

# The loop for REML-MVN 
set.seed(42)
for (i in c(1:p)){ 
  reps <- mvrnorm(n=MVNsample, theta_array[,,i], V_array[,,i])
  for (j in c(1:MVNsample)){
    M <- matrix(NA, nrow = n, ncol = n)
    lowerTriangle(M, diag = TRUE, byrow = TRUE) <- reps[j,]
    upperTriangle(M) = lowerTriangle(M, byrow = TRUE)
    AsycovM_array[,,j,i] <- M
  }
}

# Import the 1000 randomised datasets an store their M in an array to estimate nulls
# read in the data and examine the distribution of the traces of M
null_Mfilenames <- list.files(path = rando_un_cov_dir, pattern = "covpars", full.names = TRUE)

# create a function that takes in the covaraince estimates and outputs a matrix
vec_est_2mat <- function(x){
  x <- c(x)
  m <- matrix(NA, nrow=n, ncol=n)
  lowerTriangle(m, diag=TRUE, byrow = TRUE) <- x
  upperTriangle(m) = lowerTriangle(m, byrow = TRUE)
  return(m)
}

# use a combine function called acomb that uses the abind function 
# from the abind package to combine the matrices generated by the cluster

# Create a function to read in data frames and convert data into a matrix, saved in an array
rando_M_import <- function(filenum){
  # Select the among-line variance component, "Line"
  null_M_df <- fread(null_Mfilenames[filenum])[Subject == "Line"]
  rand_rep_vec <- sort(unique(null_M_df$RandRep))
  ngens = p/2
  null_array  <- array(NA, dim=c(n, n, length(rand_rep_vec), p))
  dimnames(null_array) <- list(traitnumber, traitnumber,rand_rep_vec, names(M_list))
  for (j in rand_rep_vec) {
    num_name_rep <- as.character(j)
    # Import M matrices for the small treatment (p:1 to 6)
    treat_s_M <- lapply(seq_along(1:ngens), function(x){
      vec_est_2mat(null_M_df[RandRep == j & Treat == 1 & Gen == x]$Estimate)
    })
    
    null_array[,,num_name_rep,1:ngens] <- array(unlist(treat_s_M), 
                                     dim = c(n, n, ngens))
    
    # Import M matrices for the large treatment (p:7 to 12)
    treat_b_M <- lapply(seq_along(1:ngens), function(x){
      vec_est_2mat(null_M_df[RandRep == j & Treat == 3 & Gen == x]$Estimate)
    })
    null_array[,,num_name_rep,c(ngens+1):p] <- array(unlist(treat_s_M), 
                                                     dim = c(n, n, ngens))
  }
  return(null_array)
}

n_cores <- detectCores()-1
cluster <- makeCluster(n_cores)
doParallel::registerDoParallel(cluster)
acomb <- function(...) abind(..., along=3)

# null_M_array <- 
null_M_array <- foreach (i = 1:length(null_Mfilenames), 
                         .packages = c('data.table', 'foreach','gdata'), 
                         .combine='acomb', .multicombine=TRUE) %dopar% {
                           rdm_null_array <- rando_M_import(i)
                           rdm_null_array
                         }

parallel::stopCluster(cluster)

dim(null_M_array)

# general range for 90% CI functions
rangeFunc90 <- function(x){
  n <- length(x)
  lo <- as.numeric(quantile(x, c(.05)))
  hi <- as.numeric(quantile(x, c(0.95)))
  return(c(n, lo, hi))
}

# general range for 90% CI functions
rangeFunc95 <- function(x){
  n <- length(x)
  lo <- as.numeric(quantile(x, c(.025)))
  hi <- as.numeric(quantile(x, c(0.975)))
  return(c(n, lo, hi))
}


# Section 2. M matrices eigenanalyses ------------------------------------------------------
# Output the M matrices and eigenanalyses tables with 90% CI for supplementary
# Create a figure contrasting M matrices across the 12 populations against univariaste estimates

# Generate the coordinates of the 6 x 6 matrix
comb <- data.frame(NULL)
l = 0; e = c(1:n)
while (l < 6) {
  comb <- rbind(comb, cbind(rep((l+1),(n-l)), e[(l+1):n]))
  l = l + 1
} 

# Calculated the CI for each element in M using the 10,000 REML-MVN AsycovM_array
M_tab<-data.frame(NULL)
for (mat in 1:p) {
  for (j in 1:nrow(comb)) {
    vl<-sprintf("%.3f",M_list[[mat]][comb[j,1],comb[j,2]])
    vl_ci<-paste(sprintf("%.3f",rangeFunc90(AsycovM_array[comb[j,1],comb[j,2],,mat])[2:3]), collapse="; ")
    M_tab<-rbind(M_tab,cbind(names(M_list)[[mat]],comb[j,1],comb[j,2], vl,vl_ci))
  }
}
setDT(M_tab)
M_tab[, c("Treat","Gen") := data.table(str_split_fixed(V1,"_",2))]
M_tab <- melt(M_tab[,-1],id.var=c("Treat", "Gen", "V2", "V3"), measure.vars = c("vl","vl_ci") )
dcast(M_tab, Gen+variable+V3~Treat+V2, value.var=c("value"))->M_tab
setorderv(M_tab, c("Gen", "V3"), c(1,1))
M_tab[is.na(M_tab)]<-""


# write.table(M_tab, file =paste(".",outdir_tab,"M_tab_CI.txt", sep="/"),
#           row.names = FALSE, quote = TRUE, sep="\t")

TraceM<-data.frame(NULL)
for (i in 1:p) {
  TraceM<-rbind(TraceM, tr(M_list[[i]]))
}

n_cores <- detectCores()-1
cluster <- makeCluster(n_cores)
doParallel::registerDoParallel(cluster)
acomb <- function(...) abind(..., along=4)

# Get CIs for Eigenvectors
stattime =Sys.time()
tmpVec <- lapply(M_list, function(x) t(eigen(x)$vector))
EigenMCI_array <- foreach (i = 1:MVNsample,  
                           .combine='acomb', .multicombine=TRUE) %dopar% {
                             tmpAry <- sapply(seq_along(1:p), function(x){
                               diag(tmpVec[[x]] %*% AsycovM_array[,,i,x] %*% t(tmpVec[[x]]))
                             })
                             array(c(tmpAry), dim=c(1,n,p))
                           }
end=Sys.time()-stattime; end #Time difference of 13.30392 secs
parallel::stopCluster(cluster)
dim(EigenMCI_array)

EigenM_CI_tab<-data.frame(NULL)
for (mat in 1:p) {
  for (i in 1:n) {
    EigenM_CI_tab <- rbind(EigenM_CI_tab, 
                           c(names(M_list)[[mat]],paste0("e",i),
                             eigen(M_list[[mat]])$values[i], 
                             rangeFunc90(EigenMCI_array[,i,mat,])[2:3]))
  }
}
colnames(EigenM_CI_tab)<-c("M", "eigV", "Lambda", "Lo", "Hi")
EigenM_CI_tab<-cbind(EigenM_CI_tab[,1:2],
      apply(EigenM_CI_tab[,3:5], 2, function(x) sprintf("%.3f", as.numeric(x))))
setDT(EigenM_CI_tab)
EigenM_CI_tab[, CI:= paste(Lo, Hi, sep = "; "), by=c("M", "eigV")]
EigenM_CI_tab[, c("Treat","Gen") := data.table(str_split_fixed(M,"_",2))]
EigenM_CI_tab->EigenM_CI_tab_plot

EigenM_CI_tab<-melt(EigenM_CI_tab[,.(Treat, Gen, eigV, Lambda, CI)],
                    id.var=c("Treat", "Gen", "eigV"), measure.vars = c("Lambda","CI"))
dcast(EigenM_CI_tab, Treat+Gen+variable~eigV, value.var=c("value"))->EigenM_CI_tab
EigenM_CI_tab[,Gen2:=fcase(variable=="Lambda", Gen,
                           default = "")]
EigenM_CI_tab[, Treat2:=fcase(Treat==1 & Gen==1 & variable=="Lambda", "Small",
                              Treat==3 & Gen==1 & variable=="Lambda", "Large",
                              default = "")]

# write.table(EigenM_CI_tab[,c(1:3,10:11,4:9)], file =paste(".",outdir_tab,"EigenM_CI_tab.txt", sep="/"),
# row.names = FALSE, quote = TRUE, sep="\t")

EigVec_M_tab<-data.frame(NULL)
for(i in 1:p){
  tmpM<-apply(eigen(M_list[[i]])$vectors, 2, function(x) sprintf("%.2f",x))
  tmpM<-cbind(rep(names(M_list)[[i]],6),c("CS","1.2","1.5","2.5","2.8","3.7"),tmpM)
  EigVec_M_tab<-rbind(EigVec_M_tab,tmpM)
  rm(tmpM)
}

as.data.frame(str_split_fixed(EigVec_M_tab$V1, "_",2))->EigVec_M_tab_head
setDT(EigVec_M_tab_head)
EigVec_M_tab_head[,Treat:=fcase(V1=="1" & V2=="1", "Small",
                                V1=="3" & V2=="1", "Large",
                                default="")]
colnames(EigVec_M_tab_head)[2]<-"Gen"
colnames(EigVec_M_tab)[2:8]<-c("Trait", paste0("e", 1:6))
EigVec_M_tab <- cbind(EigVec_M_tab_head[, .(Treat, Gen)], EigVec_M_tab[2:8])
rm(EigVec_M_tab_head)

# write.table(EigVec_M_tab, file =paste(".",outdir_tab,"EigVec_M_tab.txt", sep="/"),
# row.names = FALSE, quote = TRUE, sep="\t")

# fix the trace image
M_traces<-as.data.frame(cbind(1:12,unlist(lapply(M_list, tr))))
M_traces<-cbind(M_traces,as.data.frame(t(apply(apply(AsycovM_array, c(3,4), tr),2,rangeFunc90)))[,c(2:3)])
colnames(M_traces)<-c("pop", "trace","lo", "hi")
M_traces$x<-c(c(1:6)-0.15, c(7:12)-5.85)
setDT(M_traces)
# get get eigenvalues of M as proportion of trace
M_VarNeigenVec<-matrix(unlist(lapply(M_list, function(x) eigen(x)$value)), ncol=6, byrow = TRUE)
M_VarNeigenVec<-as.data.frame(cbind(1:12,apply(M_VarNeigenVec,2,"/", t(unlist(lapply(M_list, tr))))*100))
colnames(M_VarNeigenVec)<-(c("pop",paste0("e",1:6)))
#then variances (make sure to order them from largest to smallest!!)
M_Var<-matrix(unlist(lapply(M_list, diag)), ncol=6, byrow = TRUE)
M_Var<-t(apply(M_Var,1,function (x) sort(x, decreasing = TRUE)))
M_Var<-as.data.frame(cbind(1:12,apply(M_Var,2,"/", t(unlist(lapply(M_list, tr))))*100))
colnames(M_Var)<-c("pop",paste0("Var",1:6))
M_VarNeigenVec<-merge(M_VarNeigenVec, M_Var, by="pop")
setDT(M_VarNeigenVec)
M_VarNeigenVec<-melt(M_VarNeigenVec, id="pop", measure=patterns("e", "Var"))
colnames(M_VarNeigenVec)[2:4]<-c("VecNum", "Eig", "Var")
#plot
# svglite(filename = "C:/Users/carac/Dropbox/Analysis/Multivariate_MR/Multi_Latex_Documents/PNGs and EPS/M_traceNeigs.svg",
#     width = 7, height = 3.7, pointsize = 11,
#     bg = "white", system_fonts = "Times New Roman")

# grDevices::cairo_ps(filename =
#                       "C:/Users/carac/Dropbox/Analysis/Multivariate_MR/Multi_Latex_Documents/PNGs and EPS/M_traceNeigs.eps", 
#                     family = "Times", bg="transparent",pointsize=11, width = 7, height = 3.7)

par(mfrow=c(1,2), mar=c(4, 4.5,3, 1.5))
plot(NULL, ylim=c(0,1.4), xlim=c(0.5,6.5),bty="L",yaxt="n",
     xlab="Generations", #family = "Times New Roman",
     ylab=expression("Traces of "~bold(M)~"\u00B1 90% CI"))
axis(side=2, at=seq(0,1.4,0.2),family = "Times New Roman", labels = sprintf("%.1f",seq(0,1.4,0.2)), las=2)
M_traces[,segments(x,lo,x,hi)]
M_traces[,segments(x-0.05,lo,x+0.05,lo)]
M_traces[,segments(x-0.05,hi,x+0.05,hi)]
M_traces[, points(x,trace, pch=ifelse(pop<=6, 16,21),
                        bg=ifelse(pop<=6, "NA","white"))]
mtext("A", 3,family = "Times New Roman", outer=FALSE, cex=1.5,adj=-0.35, line=1)

plot(NULL, xlim=c(0.5,6.5), xlab="",yaxt="n",bty="L",
     ylab="",xaxt="n",yaxs="i", ylim=c(0,100))
M_VarNeigenVec[,boxplot(Eig~VecNum,col="white", outcex = 0.75, at=seq(0.85,5.85, 1), whisklty = 1,
                        yaxt="n",pch=4,xaxt="n",frame=F,boxwex=0.2, add = TRUE)]
M_VarNeigenVec[,boxplot(Var~VecNum,col="darkgrey", outcex = 0.75, at=seq(1.15,6.15, 1), whisklty = 1,
                        yaxt="n",pch=4,xaxt="n",frame=F,boxwex=0.2, add = TRUE)]
axis(side=1,at=3.5, family = "Times New Roman", "Eigenvectors", line=2, tick = FALSE)
axis(side=2,at=50, family = "Times New Roman", "Proportion of among-line variance", line=2, tick = FALSE)
axis(side=2, family = "Times New Roman",at=seq(0,100, 20),paste0(seq(0,100, 20),"%"), las=2)
axis(side=1,family = "Times New Roman",at=1:6, as.expression(lapply(1:6, function(i)bquote(italic("e")[.(i)]))))
mtext("B", 3, family = "Times New Roman", outer=FALSE, cex=1.5,adj=-0.35, line=1)
# dev.off()

# Section 3. Residual Varaince (not included in MS) -----------------
# Get observed R matrices from REML
# create a function to import M and convert the column list into a 6 x 6 matrix
import_R <- function(mlist_name){
  Rdf <- fread(paste(un_cov_dir,mlist_name, sep = "/"))[Subject == "Anim(Trea*Line*Vial)",.(Estimate)] 
  R <- matrix(NA, nrow = n , ncol = n)
  upperTriangle(R, diag = TRUE) <- Rdf$Estimate
  lowerTriangle(R) = upperTriangle(R, byrow = TRUE)
  return(R)
}
R_list <- lapply(Mfilenames, import_R)

names(R_list) <- gsub(".csv","", gsub("Covparms","", Mfilenames)) 


VR_list <- lapply(Vfilenames, function (x) 
  as.matrix(fread(paste(un_cov_dir, x, sep = "/"))[c(43:63), c(45:65)]))

names(VR_list) <- gsub(".csv","", gsub("Asycov","", Vfilenames)) 


#create the R and (R) V arrays and theta (R) arrays
R_array <- array(NA, dim = c(n, n, p))
for (i in 1:p){R_array[,,i] <- R_list[[i]]}
dimnames(R_array) <- list(traitnumber, traitnumber, names(M_list))

# Get asycov for residuals
VR_array<- array(0,dim=c((n*(n+1)/2),(n*(n+1)/2),p))
for (i in 1:p){VR_array[,,i] <- VR_list[[i]]}
dimnames(VR_array)[[3]]<-names(R_list)

theta_R_list <- sapply(Mfilenames, function (x)
  fread(paste(un_cov_dir, x, sep = "/"))[Subject == "Anim(Trea*Line*Vial)",.(Estimate)])
names(theta_R_list) <- gsub(".csv", "", gsub("Covparms", "", Mfilenames))  

theta_R_array <- array(NA, dim=c((n*(n+1)/2), 1, p))
for (i in c(1:p)){theta_R_array[,,i] <- theta_R_list[[i]]}
dimnames(theta_R_array)[[3]] <- names(theta_R_list)




# Create array for 10,000 REML-MVN data
AsycovR_array<-array(0, dim=c(n, n, MVNsample, p)) 

# The loop for REML-MVN 
set.seed(42)
for (i in c(1:p)){ 
  reps_R<-mvrnorm(n=MVNsample, theta_R_array[,,i], VR_array[,,i])
  for (j in c(1:MVNsample)){
    # R
    R<-matrix(0, nrow=n, ncol=n)
    lowerTriangle(R, diag=TRUE, byrow=TRUE)<-reps_R[j,]
    upperTriangle(R) = lowerTriangle(R, byrow=TRUE)
    AsycovR_array[,,j,i]<-R
  }
}


R_tab<-data.frame(NULL)
for (mat in 1:p) {
  for (j in 1:nrow(comb)) {
    resid<-sprintf("%.3f",R_list[[mat]][comb[j,1],comb[j,2]])
    resid_ci<-paste(sprintf("%.3f",rangeFunc90(AsycovR_array[comb[j,1],comb[j,2],,mat])[2:3]), collapse="; ")
    R_tab<-rbind(R_tab,cbind(names(R_list)[[mat]],comb[j,1],comb[j,2], resid,resid_ci))
  }
}
setDT(R_tab)
R_tab[, c("Treat","Gen") := data.table(str_split_fixed(V1,"_",2))]
R_tab<-melt(R_tab[,-1],id.var=c("Treat", "Gen", "V2", "V3"), measure.vars = c("resid","resid_ci") )
dcast(R_tab, Gen+variable+V3~Treat+V2, value.var=c("value"))->R_tab
setorderv(R_tab, c("Gen", "V3"), c(1,1))
R_tab[is.na(R_tab)]<-""
# write.table(R_tab, file ="C:/Users/carac/Dropbox/Analysis/Multivariate_MR/Tables/R_tab_CI.txt",
#           row.names = FALSE, quote=TRUE, sep="\t")


# Get CIs for Eigenvectors
EigenRCI_array<-array(0, dim=c(1, n, MVNsample, p))
stattime =Sys.time()
for (mat in 1:p) {
  for (i in 1:MVNsample) {
    for (vec in 1:n) {
      tmpVec<-t(eigen(R_list[[mat]])$vector[,vec])
      EigenRCI_array[,vec,i,mat] <- tmpVec %*%
        AsycovR_array[,,i,mat] %*%
        t(tmpVec)
    }
  }
};end=Sys.time()-stattime; end #Time difference of 3.135797 mins

EigenR_CI_tab<-data.frame(NULL)
for (mat in 1:p) {
  for (i in 1:n) {
    EigenR_CI_tab <- rbind(EigenR_CI_tab,
                           c(names(R_list)[[mat]],paste0("e",i),
                             eigen(R_list[[mat]])$values[i],
                             rangeFunc90(EigenRCI_array[,i,,mat])[2:3]))
  }
}
colnames(EigenR_CI_tab)<-c("R", "eigV", "Lambda", "Lo", "Hi")
EigenR_CI_tab<-cbind(EigenR_CI_tab[,1:2],
                     apply(EigenR_CI_tab[,3:5], 2, function(x) sprintf("%.3f", as.numeric(x))))
setDT(EigenR_CI_tab)
EigenR_CI_tab[, CI:= paste(Lo, Hi, sep = "; "), by=c("R", "eigV")]
EigenR_CI_tab[, c("Treat","Gen") := data.table(str_split_fixed(R,"_",2))]
EigenR_CI_tab->EigenR_CI_tab_plot

EigenR_CI_tab<-melt(EigenR_CI_tab[,.(Treat, Gen, eigV, Lambda, CI)],
                    id.var=c("Treat", "Gen", "eigV"), measure.vars = c("Lambda","CI"))
dcast(EigenR_CI_tab, Treat+Gen+variable~eigV, value.var=c("value"))->EigenR_CI_tab
EigenR_CI_tab[,Gen2:=fcase(variable=="Lambda", Gen,
                           default = "")]
EigenR_CI_tab[, Treat2:=fcase(Treat==1 & Gen==1 & variable=="Lambda", "Small",
                              Treat==3 & Gen==1 & variable=="Lambda", "Large",
                              default = "")]

# write.table(EigenR_CI_tab[,c(1:3,10:11,4:9)], file ="C:/Users/carac/Dropbox/Analysis/Multivariate_MR/Tables/EigenR_CI_tab.txt",
# row.names = FALSE, quote = TRUE, sep="\t")

#plot CIs to confirm overlap
EigenR_CI_tab_plot[, e:=as.numeric(gsub("e","", eigV))]
as.data.table(apply(EigenR_CI_tab_plot[,c(3:5,7:9)], 2, function(x) as.numeric(x)))->EigenR_CI_tab_plot
setDT(EigenR_CI_tab_plot)
EigenR_CI_tab_plot[, ID := .I]

par(mfrow=c(1,1),mar=c(4,8,0.2,4))
plot(NULL, xlim=c(0, 2),yaxt="n",ylab="", xlab="Variance", ylim=c(1,72), bty="L")
axis(side=2,at=c(18,54), tick = FALSE,las=2, labels = c("Small", "Large"), line=3)
axis(side=2, at=seq(0,66, length.out=12), line=1.5, tick=FALSE,labels = rep(1:6,2), las=2)
axis(side=2, at=1:72, labels = rep(paste0("e", 1:6),12), cex.axis=0.5, las=2)
abline(v=0,lty=2)
for (trt in c(1,3)) {
  for (g in 1:6) {
    for (ev in 1:6){
      EigenR_CI_tab_plot[Treat==trt& Gen==g & e==ev] %T>% 
        with(segments(Lo, ID, Hi, ID, lwd=ifelse(Lo>0,2,1), 
                      col=ifelse(Lo>0 & Lambda>0, "blue",
                                 ifelse(trt==1, 1, 2)))) %>% 
        with(points(Lambda, ID,col=ifelse(Lo>0 & Lambda>0, "blue",
                                          ifelse(trt==1, 1, 2)), pch=16, cex=0.95))
    }
    
  }
}
EigenR_CI_tab_plot[, Totvar:=tr(R_list[[grep(paste(Treat,Gen,sep="_"),M_names)]]), by="ID"]
EigenR_CI_tab_plot[,PropVar:= round((Lambda/Totvar)*100,2), by="ID"]

par(mfrow=c(1,1))
EigenR_CI_tab_plot[e==1, .(PropVar, e)] %>%
  with(boxplot(PropVar ~ e,
               xlim=c(0,7), ylim=c(0,50),
               xlab="Eigenvectors",bty="L",
               ylab="Propotion of varaince (%)"))
for (i in 2:n) {
  EigenR_CI_tab_plot[e==i, .(PropVar, e)] %>% 
    with(boxplot(PropVar~e,at=i, add=TRUE,  boxwex=0.55))
}
mean(EigenR_CI_tab_plot[e==1, PropVar])
mean(EigenR_CI_tab_plot[e==2, PropVar])

# Chapter Figure
R_MVN_trace_CI<-as.data.frame(t(apply(apply(AsycovR_array[,,,], c(3,4), 
                                            function (x) tr(x)), 2, rangeFunc90)))
setDT(R_MVN_trace_CI)
colnames(R_MVN_trace_CI)<-c("n","lo","hi")
R_MVN_trace_CI[, p:=.I]
R_MVN_trace_CI$Trace<-unlist(lapply(R_list, tr))
R_MVN_trace_CI[, x:=fcase(p <=6, p-0.15,
                          p>6, p - 5.85, default=0)]

R_eigenVectors<-matrix(unlist(lapply(R_list, 
                                     function (x) eigen(x)$value)), ncol=6, byrow=TRUE)
R_eigenVectors<-as.data.frame(apply(R_eigenVectors,2,"/", t(unlist(lapply(R_list, tr))))*100)
setDT(R_eigenVectors)
colnames(R_eigenVectors)<-paste0("e", 1:6)
R_eigenVectors[,p:=.I]
melt(R_eigenVectors, id.vars = "p")->R_eigenVectors
R_eigenVectors$variable<-factor(R_eigenVectors$variable, levels=paste0("e", 1:6))

# svglite(filename = "C:/Users/carac/Dropbox/Analysis/Multivariate_MR/Multi_Latex_Documents/PNGs and EPS/R_traceNeigs.svg",
#     width = 7, height = 3.6, pointsize = 11,
#     bg = "white", system_fonts = "Times New Roman")

grDevices::cairo_ps(filename =
                      "C:/Users/carac/Dropbox/Analysis/Multivariate_MR/Multi_Latex_Documents/PNGs and EPS/R_traceNeigs.eps", 
                    family = "Times", bg="transparent",pointsize=11, width = 7, height = 3.6)
# R traces
par(mfrow=c(1,2), mar=c(4, 4.5,3, 1.5))
plot(NULL, ylim=c(0,6), xlim=c(0.5,6.5),bty="L",yaxt="n",
     xlab="Generations", family = "Times New Roman",
     ylab=expression("Traces of "~bold(R)~"\u00B1 90% CI"))
axis(side=2, at=seq(0,6),family = "Times New Roman", labels = sprintf("%.1f",seq(0,6)), las=2)
R_MVN_trace_CI[,segments(x,lo,x,hi)]
R_MVN_trace_CI[,segments(x-0.05,lo,x+0.05,lo)]
R_MVN_trace_CI[,segments(x-0.05,hi,x+0.05,hi)]
R_MVN_trace_CI[, points(x,Trace, pch=ifelse(p<=6, 16,21),
                        bg=ifelse(p<=6, "NA","white"))]
mtext("A", 3,family = "Times New Roman", outer=FALSE, cex=1.5,adj=-0.35, line=1)

R_eigenVectors %>% 
  with(plot(NULL, xlim=c(0.5,6.5), 
            xlab="",yaxt="n",bty="L",
            ylab="",xaxt="n", 
            yaxs="i", ylim=c(0,40)))
R_eigenVectors[,boxplot(value~variable,col="white",
                        yaxt="n",pch=4,xaxt="n",frame=F,boxwex=0.35, add = TRUE)]
axis(side=1,at=3.5, family = "Times New Roman", "Eigenvectors", line=2, tick = FALSE)
axis(side=2,at=20, family = "Times New Roman", "Proportion of residual variance", line=2, tick = FALSE)
axis(side=2, family = "Times New Roman",at=seq(0,40, 5),paste0(seq(0,40, 5),"%"), las=2)
axis(side=1,family = "Times New Roman",at=1:6, as.expression(lapply(1:6, function(i)bquote(italic("e")[.(i)]))))
mtext("B", 3, family = "Times New Roman", outer=FALSE, cex=1.5,adj=-0.35, line=1)
# dev.off()



# Section 4. Krzanowksi's Common Subspaces, H ------------------
# Part A. Implement Krzanowksi's method, and generate eigenanalysis and null distributions 
# First, get 10,000 REML-MVN h matrices

# using the first 3 largest eigenvectors from M to determine in H
# then saving the eigenvectors of H
H_vecs <- evolqg::KrzSubspace(M_list, 3)$k_eVec_H

n_cores <- detectCores()-1
cluster <- makeCluster(n_cores)
doParallel::registerDoParallel(cluster)
acomb <- function(...) abind(..., along=3)

remlmvn_H_array <- foreach (i = 1:MVNsample, 
                         .packages = c('evolqg'), 
                         .combine='acomb', .multicombine=TRUE) %dopar% {
                           temp_rml_list <- asplit(AsycovM_array[,,i,], 3)
                           evolqg::KrzSubspace(temp_rml_list, 3)$H
                         }

null_H_array <- foreach(i = 1:nullnumber, 
                        .packages = c('evolqg'), 
                        .combine='acomb', .multicombine=TRUE) %dopar% {
                          temp_null_list <- asplit(null_M_array[,,i,], 3)
                          evolqg::KrzSubspace(temp_null_list, 3)$H
                        }

parallel::stopCluster(cluster)

dim(remlmvn_H_array)
dim(null_H_array)

# Get CIs for H eigenvalues
# Calculate the eigenvalues for each of the 10,000 H matrices for the CIs
remlmvn_H_eigvals <- foreach(i = 1:MVNsample, 
                             .combine='rbind') %do% {
                               eigen(remlmvn_H_array[,,i])$value[1:3]
                             }
H_eigval_CI <- t(apply(remlmvn_H_eigvals, 2, function(x) rangeFunc95(x)[2:3]))

H_eigval_CI <- cbind(1:3, evolqg::KrzSubspace(M_list, 3)$k_eVals_H,
                     H_eigval_CI)
colnames(H_eigval_CI)<- c("vec_num","H_eigval", "CI_lo", "CI_up")

H_eigval_CI <- setDT(as.data.frame(H_eigval_CI))

H_eigval_CI$Est <- "Obs"

# Get mean and CIs for the null
avg_null_H <- apply(null_H_array, 1:2, mean)
avg_null_H_vecs <- eigen(avg_null_H)$vectors[,1:3]

null_H_eigvals <- laply(asplit(null_H_array, 3), 
                           function(mat) diag(t(avg_null_H_vecs) %*% mat %*% avg_null_H_vecs))
null_H_CI <- t(apply(null_H_eigvals, 2, function(x) rangeFunc95(x)[2:3]))

null_H_CI <- cbind(as.integer(row.names(null_H_CI)),
                   eigen(avg_null_H)$values[1:3], null_H_CI)
colnames(null_H_CI)<- c("vec_num","H_eigval", "CI_lo", "CI_up")

null_H_CI <- setDT(as.data.frame(null_H_CI))

null_H_CI$Est <- "Null"

# Recreate Dave's Figure 3.
par(mfrow=c(1,1), mar=c(4, 4.5,3, 1.5))
plot(NULL, ylim=c(0,12), xlim=c(0.5,3.5),bty="L",xaxt="n",las=2,
     xlab=expression("Eigenvectors of"~bold(H)), #family = "Times New Roman"
     ylab=expression("Eigenvalues of"~bold(H)~"\u00B1 95% CI"))
H_eigval_CI[,segments(vec_num-0.25,CI_lo,vec_num-0.25,CI_up)]
H_eigval_CI[,segments(vec_num-0.20,CI_lo,vec_num-0.30,CI_lo)]
H_eigval_CI[,segments(vec_num-0.20,CI_up,vec_num-0.30,CI_up)]
H_eigval_CI[,points(vec_num-0.25, H_eigval, pch=16)]
axis(side=1,family = "Times New Roman",at=1:3, as.expression(lapply(1:3, function(i)bquote(italic("h")[.(i)]))))
null_H_CI[,segments(vec_num+0.25,CI_lo,vec_num+0.25,CI_up, lty=2)]
null_H_CI[,segments(vec_num+0.20,CI_lo,vec_num+0.30,CI_lo, lty=2)]
null_H_CI[,segments(vec_num+0.20,CI_up,vec_num+0.30,CI_up, lty=2)]
null_H_CI[,points(vec_num+0.25, H_eigval)]
legend("bottomright", c("Observed", "Randomised"), bty="n", 
       pch=c(16,21), lty=c(1,2))


# Create a table
H_tab <- rbind(H_eigval_CI,null_H_CI)

H_tab <- cbind(H_tab$Est, paste0("h",H_tab$vec_num),
               apply(H_tab[,2:4], 2, function(x) sprintf("%.2f", as.numeric(x))))

colnames(H_tab) <- c("Est", "eigV", "Lambda", "Lo", "Hi")
H_tab <- as.data.table(H_tab)
H_tab[, CI:= paste(Lo, Hi, sep = "; "), by=c("Est", "eigV")]

H_tab <- melt(H_tab[,.(Est, eigV, Lambda, CI)], 
              measure.vars = c("Lambda","CI"))

H_tab <- dcast(H_tab,Est+variable~eigV)

# write.table(H_tab, file =paste(".",outdir_tab,"EigH_tab.txt", sep="/"),
# row.names = FALSE, quote = TRUE, sep="\t")


# Part B. Explore the amount of among-line variance each h eigenvector captures
# Looking at Krzanowksi's H eigvec projections through M 
Var_HeigVecs<-array(NA, dim=c(1, p ,3))
for(j in c(1:3)){#number of eigenvectors of H
  for (i in c(1:p)) {
    Var_HeigVecs[,i,j] <- t(H_vecs[,j]) %*% M_array[,,i] %*% t(t(H_vecs[,j]))
  }
}

# Get CIS
AsycovHeigVecs_array<-array(NA, dim=c(1,1, MVNsample, p, 3))
for(j in c(1:3)){#number of eigenvectors of H
  for (i in c(1:p)) {
    for (k in 1:MVNsample){
      AsycovHeigVecs_array[,,k,i,j] <- t(H_vecs[,j]) %*% AsycovM_array[,,k,i] %*% t(t(H_vecs[,j]))
    }
  }
}



# setEPS()
# postscript("C:/Users/carac/Dropbox/Analysis/Multivariate_MR/Multi_Latex_Documents/PNGs and EPS/HvecsThruMnR.eps", 
#            family = "Times", pointsize=14, width =10, height = 3.25)


# svglite(filename = "C:/Users/carac/Dropbox/Analysis/Multivariate_MR/Multi_Latex_Documents/PNGs and EPS/HvecsThruMnR.svg",
#     width = 6, height = 5.9, pointsize = 10,
#     bg = "white", system_fonts = "Times New Roman")

# grDevices::cairo_ps(filename =
#                       "C:/Users/carac/Dropbox/Analysis/Multivariate_MR/Multi_Latex_Documents/PNGs and EPS/HvecsThruMnR.eps",
#                     family = "Times", bg="transparent",pointsize=10, width = 6, height = 5.9)

as.data.frame(cbind(rep(1:12, 6), c(rep("M", 36),rep("R", 36)),
      rep(rep(1:3,each=6),2),rep(c(seq(1.25,11.25,2), seq(1.75,11.75,2)),6)))->H_thru_MnR
setDT(H_thru_MnR)
H_thru_MnR[,V1:=as.numeric(V1)][,V3:=as.numeric(V3)][,V4:=as.numeric(V4)]
colnames(H_thru_MnR)<-c("p", "parm","Hvec","x")
H_thru_MnR[, pop:=c(unlist(names(M_list)[p])), by=.I]
H_thru_MnR[,c("var","Lo", "Hi"):=0]

for(j in c(1:3)){#number of eigenvecotrs of H
  for (i in c(1:p)) {
    H_thru_MnR[p==i & Hvec==j & parm=="M", var:= Var_HeigVecs[,i,j]]
    tempCI=rangeFunc90(AsycovHeigVecs_array[,,,i,j])
    H_thru_MnR[p==i & Hvec==j & parm=="M", `:=`(Lo=tempCI[2] , Hi=tempCI[3])]
    H_thru_MnR[p==i & Hvec==j & parm=="R", var:= Var_HeigVecs_R[,i,j]]
    tempCI=rangeFunc90(AsycovHeigVecs_array_R[,,,i,j])
    H_thru_MnR[p==i & Hvec==j & parm=="R", `:=`(Lo=tempCI[2] , Hi=tempCI[3])]
  }
};rm(tempCI)
H_thru_MnR[,pchy:=fcase(p %in% 1:6, 16,
                        default= 21)]

par(mfrow=c(3,3),mar = c(4, 5, 4, 0.5))
for (prm in c("M", "R")) {
  for(j in c(1:3)){
    if (prm=="M"){
      plot(NULL, xlim=c(0.5,12.5), ylim=c(-0.1, 0.75), bty="L",xaxt="n",
           ylab=ifelse(j==1, "Among-line variance \u00B1 90% CI", ""),
           xlab="Generations", yaxt="n", family = "Times New Roman", 
           main=bquote(bolditalic("h")[.(j)]))
      abline(h=0, lty=2)
      axis(side=1, at=seq(1.5,11.5,2), labels = 1:6, family = "Times New Roman")
      axis(side=2, at=seq(-0.1,0.7,0.1), labels = seq(-0.1,0.7,0.1), las=2, family = "Times New Roman")
      mtext(LETTERS[j], side=3, adj=-0.2, line=1.3,cex=1.2, family = "Times New Roman")
    }else{
      plot(NULL, xlim=c(0.5,12.5), ylim=c(0.4, 1.8), bty="L",xaxt="n",
           ylab=ifelse(j==1, "Residual variance \u00B1 90% CI", ""),
           xlab="Generations", yaxt="n", family = "Times New Roman", 
           main=bquote(bolditalic("h")[.(j)]))
      axis(side=1, at=seq(1.5,11.5,2), labels = 1:6, family = "Times New Roman")
      axis(side=2, at=seq(0.4,1.8,0.2), labels = sprintf("%1.1f",seq(0.4,1.8,0.2)), las=2, family = "Times New Roman")
      mtext(LETTERS[j+3], side=3, adj=-0.2, line=1.3,cex=1.2, family = "Times New Roman")
    }
    for(i in 1:p){
    H_thru_MnR[Hvec==j & parm==prm & p==i] %T>%
        with(segments(x,Lo,x, Hi)) %T>% 
        with(points(x,var, pch=pchy, bg=ifelse(pchy==21, "white", NA))) %T>% 
        with(segments((x-0.1), Lo,(x+0.1), Lo)) %>% 
        with(segments((x-0.1), Hi,(x+0.1), Hi))
    }
  }
}
H_thru_MnR_w<- H_thru_MnR
H_thru_MnR_w[, c("Treat","Gen") := data.table(str_split_fixed(pop,"_",2))]
H_thru_MnR_w[,c("x", "pop", "Treat"):= NULL]
H_thru_MnR_w[,Gen:=as.numeric(Gen)]
dcast(H_thru_MnR_w, p+Gen+Hvec+pchy~parm, value.var = c("var","Lo","Hi"))->H_thru_MnR_w
#Spearman's or Pearsons?
# for(j in c(1:3)){
#   cat(paste0("Hvec ",j, "\n"))
#   cat(ifelse(shapiro.test(H_thru_MnR_w[Hvec==j]$var_M)$p.value > 0.05, "M normal", "M not normal"), "\n")
#   cat(ifelse(shapiro.test(H_thru_MnR_w[Hvec==j]$var_R)$p.value > 0.05, "R normal", "R not normal"), "\n")
#   cat("____\n\n")
# }
for (j in c(1:3)) {
  plot(NULL, xlim=c(0.6,1.5), ylim=c(0, 0.5), bty="L",xaxt="n",
       ylab=ifelse(j==1, "Among-line variance", ""),
       xlab="Residual variance", yaxt="n", family = "Times New Roman", 
       main=bquote(bolditalic("h")[.(j)]))
  axis(side=1, at=seq(0.6,1.5,0.2), family = "Times New Roman", labels = sprintf("%1.1f",seq(0.6,1.5,0.2)))
  axis(side=2, at=seq(0,0.5,0.1), family = "Times New Roman", labels = seq(0,0.5,0.1), las=2)
  mtext(LETTERS[j+6], side=3, adj=-0.2, line=1.3,cex=1.2, family = "Times New Roman")
  for(i in 1:p){
    H_thru_MnR_w[Hvec==j & p==i] %T>% 
      with(points(var_R,var_M, pch=pchy, bg=ifelse(pchy==21, "white", NA))) %>% 
      with(text(var_R,var_M,Gen,pos=2,cex=0.95, family = "Times New Roman"))
  }
  if (j!=4){#3
  tempCorvals<-H_thru_MnR_w[Hvec==j, cor.test(var_M,var_R, method="pearson")][c(4,3)]
  mtext(bquote("(Pearson's "~italic(r)~"="~.(sprintf("%1.2f",tempCorvals$estimate))~
                 ","~italic(P)~"="~.(sprintf("%1.2f",tempCorvals$p.value))~" )"),
        line=-0.1, family = "Times New Roman", cex=0.65)
  }else{
      tempCorvals<-H_thru_MnR_w[Hvec==j, cor.test(var_M,var_R, method="spearman")][c(4,3)]
  mtext(bquote("(Spearman's "~italic(rho)~"="~.(sprintf("%1.2f",tempCorvals$estimate))~
                 ","~italic(P)~"="~.(sprintf("%1.2f",tempCorvals$p.value))~" )"),
        line=-0.1, family = "Times New Roman", cex=0.65)
  }
       
}
 # dev.off()


# For chapter - get H
# apply(evolqg::KrzSubspace(R_list, 3)$H, 2, function(x) sprintf("%.3f", x)) %>% View()
```
```{r Eigentensor Analysis}
# Tensor Analysis for REML estimates ONLY ---------------------------
# Construction of M covariance tensor S using REML estimates [Dave's Method]
neigten<-n*(n+1)/2
REML_S<-array(NA, dim=c(neigten,neigten))
dimnames(REML_S)<-list(paste0("e", 1:neigten), paste0("e", 1:neigten))
REML_M<-M_array
REML_Mvarmat<-t(apply(REML_M, 3, diag))
REML_Mcovmat<-t(apply(REML_M, 3, lowerTriangle))
REML_S[1:n, 1:n]<- cov(REML_Mvarmat,REML_Mvarmat) #upper left quarter of S
REML_S[(n+1):neigten, (n+1):neigten]<-2*cov(REML_Mcovmat,REML_Mcovmat)#lower right quarter of S
REML_S[1:n, (n+1):neigten]<-sqrt(2)*cov(REML_Mvarmat, REML_Mcovmat)
REML_S[(n+1):neigten, 1:n]<-sqrt(2)*cov(REML_Mcovmat, REML_Mvarmat)

# Construction of M eigentensors  
#Get the S tensor Eigenvectors and Eigenvalues 
REML_S_eigvec <- eigen(REML_S)$vectors
REML_S_eigval <- eigen(REML_S)$values

# Now need to create the Eigentensors
REML_S_eigTenmat<-array(NA, dim=c(n,n,neigten))
dimnames(REML_S_eigTenmat)<-list(traitnumber, traitnumber,paste0("E", 1:neigten))
for (i in 1:neigten){
  REML_emat<-matrix(NA,n,n)
  lowerTriangle(REML_emat)<-1/sqrt(2)*REML_S_eigvec[(n+1):neigten,i]
  REML_emat<-REML_emat+t(REML_emat)
  diag(REML_emat)<-REML_S_eigvec[1:n,i]
  REML_S_eigTenmat[,,i]<-REML_emat
}

REML_S_eigten<-data.frame(NULL)
for (i in 1:neigten){REML_S_eigten<-rbind(REML_S_eigten,REML_S_eigTenmat[,,i])}

REML_S_eigtenvecs<-matrix(nrow=n*neigten, ncol=n)
REML_S_eigtenvals<-matrix(nrow=n*neigten, ncol=1)
for (i in 1:neigten){ #Using Emma's Method
  REML_S_eigtenvecs[((i-1)*n+1):(i*n),]=t(eigen(REML_S_eigten[((i-1)*n+1):(i*n),])$vectors) # eigenvectors in rows!!!!!!
  REML_S_eigtenvals[((i-1)*n+1):(i*n),1]=eigen(REML_S_eigten[((i-1)*n+1):(i*n),])$values
}


# Let's look at the Coordinates of the M matrix ---------------------------
#This is to understand how treat x gen 'populations' vary in relation to one particular independent change in genetic variances
REML_S_eigTenmat

calMcoords<-function(eigtenslist, mlist, mnames){# 
  Coords<-data.frame()
  for (e in 1:neigten){
    eigCoords<-data.frame()
    for (i in c(1:length(mlist))){
      eigCoords[i,1]<-i
      eigCoords[i,2]<-names(M_list)[i]
      eigCoords[i,3]<-paste0("C_E", e)
      eigCoords[i,4]<-frobenius.prod(eigtenslist[,,e], M_list[[i]])  
    }
    Coords<-rbind(Coords, eigCoords)
  }
  return(Coords)
}

calMcoords(REML_S_eigTenmat, M_list, names(M_list))->cords_Meigten

traceM<-data.frame(NULL)
for (i in 1:length(M_list)){
  traceM[i,1]<-i
  traceM[i,2]<-psych::tr(t(M_list[[i]]) %*% (M_list[[i]]))
}
colnames(traceM)<-c("pop", "NormGf")

cords_Meigten %>% magrittr::set_colnames(c("pop", "m_name","cord_Eig", "Coord")) %>% 
  mutate(Csqrd=Coord^2) %>% left_join(traceM, by=c("pop")) %>% 
  mutate(perC=(Csqrd/NormGf)*100) ->cords_withpercent

#For confidence intervals see after REML-MVN sampling

# Create Dave's Fig. 4 ----------------------------------------------------
# Get 10,000 S matrices from the 10,000 REML-MVN M matrix sets 
neigten<-n*(n+1)/2
MVN_S<-array(NA, dim=c(neigten,neigten,MVNsample))
dimnames(MVN_S)<-list(paste0("e", 1:neigten), paste0("e", 1:neigten))
for (i in 1:MVNsample){
  # 1. Construction of M covariance tensor
  MVN_M<-AsycovM_array[,,i,]
  MVN_Mvarmat<-t(apply(MVN_M, 3, diag))
  MVN_Mcovmat<-t(apply(MVN_M, 3, lowerTriangle))
  MVN_S[1:n, 1:n, i]<- cov(MVN_Mvarmat,MVN_Mvarmat) #upper left quarter of S
  MVN_S[(n+1):neigten, (n+1):neigten, i]<-2*cov(MVN_Mcovmat,MVN_Mcovmat)#lower right quarter of S
  MVN_S[1:n, (n+1):neigten,i]<-sqrt(2)*cov(MVN_Mvarmat, MVN_Mcovmat)
  MVN_S[(n+1):neigten, 1:n, i]<-sqrt(2)*cov(MVN_Mcovmat, MVN_Mvarmat)
}

#Import in our best estimate (REML) of Eigentensors (E)
REML_eigten<- REML_S_eigTenmat

#Get the eigen values and the eigenvectors 
#already calculated REML_S
REML_S_val <- eigen(REML_S)$values
REML_S_vec <- eigen(REML_S)$vectors


MVN_S_val <- matrix(NA, MVNsample, neigten)
colnames(MVN_S_val) <- paste("E", 1:neigten, sep="")
for (i in 1:MVNsample){
  for(j in 1:neigten){
    MVN_S_val[i,j] <- t(REML_S_vec[,j]) %*% MVN_S[,,i] %*% REML_S_vec[,j]
  }
}

# Plot the 90% CIs for eigtensors -----------------------------------------
#Add REML-MVN CI
EigtenCIs<-as.data.frame(NULL)
for (p in 1:neigten){
  EigtenCIs<-rbind(EigtenCIs, c(p, rangeFunc90(MVN_S_val[,p])))
}
colnames(EigtenCIs)<-c("Eigten_Num", "n", "lowCI", "upCI")
EigtenCIs %<>% mutate(Eigten=paste0("E", EigtenCIs$Eigten_Num))

#add reml estimates
REML_S_eigval %>% as.data.frame %>% 
  magrittr::set_colnames("REML_Var") %>%  
  mutate(Eigten_Num=row.names(.)) %>% 
  mutate(Eigten_Num=as.numeric(Eigten_Num)) %>%
  left_join(EigtenCIs, by=c("Eigten_Num"))->plot_EigtenCIs

# Get Major axes of Eigentensor 1 e11
Project_major_EigtenVecs<-matrix(NA, nrow=0, ncol=6)
for (evec in c(1,6,7,13,18)) {
  e_vec<-REML_S_eigtenvecs[evec,]
  enum=ifelse(evec<7, (10+evec), ifelse(evec>10, evec+18, 21))
  # Get CIs for M
  e_vec_thrmat<-t(apply(apply(AsycovM_array, c(3,4), function (x) t(e_vec) %*% x %*% e_vec),2,function (t) rangeFunc90(t)))
  e_vec_thrmat<-cbind(rep(13,12),rep(enum,12), apply(M_array, 3, function (x) t(e_vec) %*% x %*% e_vec), e_vec_thrmat)
  Project_major_EigtenVecs<-rbind(Project_major_EigtenVecs,e_vec_thrmat)
  # Get CIs for R
  e_vec_thrmat<-t(apply(apply(AsycovR_array, c(3,4), function (x) t(e_vec) %*% x %*% e_vec),2,function (t) rangeFunc90(t)))
  e_vec_thrmat<-cbind(rep(18,12),rep(enum,12), apply(R_array, 3, function (x) t(e_vec) %*% x %*% e_vec), e_vec_thrmat)
  Project_major_EigtenVecs<-rbind(Project_major_EigtenVecs,e_vec_thrmat)
}; rm(e_vec_thrmat)
colnames(Project_major_EigtenVecs)<-c("Parm", "vec", "var", "n", "loCI", "upCI")
matInfo<-rownames(Project_major_EigtenVecs)
Project_major_EigtenVecs<-as.data.frame(Project_major_EigtenVecs)
setDT(Project_major_EigtenVecs)
Project_major_EigtenVecs[, Parm:=LETTERS[Parm]]
Project_major_EigtenVecs$matInfo<-matInfo
Project_major_EigtenVecs[, c("Treat","Gen") := data.table(str_split_fixed(matInfo,"_",2))]
Project_major_EigtenVecs[, x:=fcase(Treat=="1", as.numeric(Gen)-0.15,
                                    Treat=="3", as.numeric(Gen)+0.15,default=NA)]
dcast(Project_major_EigtenVecs[,.(Treat, Gen, Parm, vec, var)],
      vec+Treat+Gen~Parm, value.var="var")->Project_major_EigtenVecs_w
# # Plot the variances explained by the Eigentensors with the REML-MVN CIs
# setEPS()
# postscript("C:/Users/uqcconra/Dropbox/3_ MR data/Analysis/Univariate_MR/Latex_Documents/DaveFig_nMR.eps",
#            family = "Times", pointsize=14, width = 7.5, height = 8)

# emf("DaveFig_nMR.emf", pointsize=14, width = 7, height = 7.7)

# svglite(filename = "C:/Users/carac/Dropbox/Analysis/Multivariate_MR/Multi_Latex_Documents/PNGs and EPS/DaveFig_nMR.svg",
#     width = 6, height = 8.1, pointsize = 10,
#     bg = "white", system_fonts = "Times New Roman")
par(mfrow=c(1,1))
layout(matrix(c(rep(1,8),2,2,3,3,2,2,3,3,4,4,5,5,4,4,5,5,6,6,7,7,6,6,7,7), 8, 4, byrow = TRUE))
# layout.show(5)
Eigtenlim<-11
par(mar = c(4, 14, 1.5, 12))
plot(NULL, type="n", xlab="", yaxs="i", ylab="Among-Line Variance", 
     main="",  xaxt="n", yaxt="n",family = "Times New Roman",
     xlim=c(1,Eigtenlim), ylim=c(0,0.06), bty="L")
epsilon=0.1
for (i in c(1:Eigtenlim)){
  plot_EigtenCIs %>% filter(Eigten_Num==i)->subEig
  with(subEig, segments(Eigten_Num, lowCI, Eigten_Num, upCI))
  with(subEig, segments(Eigten_Num-epsilon, lowCI, Eigten_Num+epsilon, lowCI))
  with(subEig, segments(Eigten_Num-epsilon, upCI, Eigten_Num+epsilon, upCI))
  #with(subEig, points(Eigten_Num, Median))
  with(subEig, points(Eigten_Num, REML_Var, pch=16))
}
axis(2, at=seq(0,0.06, 0.01), labels=sprintf("%.2f", seq(0,0.06, 0.01)), las=2,family = "Times New Roman")
axis(1, at=c(1:Eigtenlim), family = "Times New Roman",
     labels=as.expression(lapply(1:Eigtenlim, function(i)bquote(bold("E")[.(i)]))))
mtext("A", 3, outer=FALSE, cex=1.25,family = "Times New Roman", adj=-0.21, line=0.1)

par(mar=c(4.5,6,3.5,1)); t=2
for (p in c("M", "R")) {
  for (evec in c(11,16)) {
    if (p=="M"){
      plot(NULL, xlim=c(0.5, 6.5),bty="L",yaxt="n",xlab="Generations",
           ylim=c(-0.05, 1), main="",
           family = "Times New Roman",
           ylab="Among-line varaince \u00B1 90% CI")
      axis(side=2, at=seq(0,1,0.2),family = "Times New Roman", labels = sprintf("%.2f",seq(0,1,0.2)), las=2)
      abline(h=0, lty=2)
      mtext(bquote(italic("e")[.(evec)]),3, adj=0.03, line = 0,family = "Times New Roman")
    } else{
      plot(NULL, xlim=c(0.5, 6.5),bty="L",yaxt="n",xlab="Generations",
           ylim=c(0.5, 2), main="",
           family = "Times New Roman",
           ylab="Residual variance \u00B1 90% CI")
      axis(side=2, at=seq(0.6,2,0.2),family = "Times New Roman", labels = sprintf("%.2f",seq(0.6,2,0.2)), las=2)
      mtext(bquote(italic("e")[.(evec)]),3, adj=0.03, line = 0,family = "Times New Roman")
    } 
    Project_major_EigtenVecs[Parm==p & vec == evec,segments(x,loCI,x, upCI), by=.I]
    Project_major_EigtenVecs[Parm==p & vec == evec,segments(x-0.05,loCI,x+0.05, loCI), by=.I]
    Project_major_EigtenVecs[Parm==p & vec == evec,segments(x-0.05,upCI,x+0.05, upCI), by=.I]
    Project_major_EigtenVecs[Parm==p & vec == evec,
                             points(x,var, pch=21, bg=ifelse(Treat=="1","black", "white")), by=.I]
    mtext(LETTERS[t], 3, outer=FALSE, cex=1.25, adj=-0.21,family = "Times New Roman", line=1.55);t=t+1
  }
}
Project_major_EigtenVecs_w[vec==11] %T>% 
  with(plot(R,M, bty="L", xlim=c(0.8, 1.8),ylim=c(0,0.5), yaxt="n",
            xlab="Residual variance", ylab="Among-line variance",
            pch=ifelse(Treat==1, 16,21),family = "Times New Roman")) %>% 
  with(text(R,M,Gen,pos=2,cex=0.95, family = "Times New Roman"))
axis(side=2, at=seq(0,0.5, 0.1),family = "Times New Roman", 
     labels = sprintf("%1.1f",seq(0,0.5, 0.1)), las=2)
mtext(bquote(italic("e")[.(11)]),3, line = 0.95,family = "Times New Roman")
tempCorvals<-Project_major_EigtenVecs_w[vec==11, cor.test(R,M, method="pearson")][c(4,3)]
mtext(bquote("(Pearson's "~italic(r)~"="~.(sprintf("%1.2f",tempCorvals$estimate))~
               ","~italic(P)~"="~.(sprintf("%1.2f",tempCorvals$p.value))~" )"),
      line=-0.1, family = "Times New Roman", cex=0.65)
mtext(LETTERS[6], 3, outer=FALSE, cex=1.25, adj=-0.21,family = "Times New Roman", line=1.55)
Project_major_EigtenVecs_w[vec==16] %T>% 
  with(plot(R,M, bty="L", xlim=c(0.7, 1.4),ylim=c(0,0.5),yaxt="n",
            xlab="Residual variance", ylab="Among-line variance",
            pch=ifelse(Treat==1, 16,21),family = "Times New Roman")) %>% 
  with(text(R,M,Gen,pos=2,cex=0.95, family = "Times New Roman"))
mtext(bquote(italic("e")[.(16)]),3, line = 0.95,family = "Times New Roman")
axis(side=2, at=seq(0,0.5, 0.1),family = "Times New Roman", 
     labels = sprintf("%1.1f",seq(0,0.5, 0.1)), las=2)
tempCorvals<-Project_major_EigtenVecs_w[vec==11, cor.test(R, M, method="spearman")][c(4,3)]
mtext(bquote("(Spearman's "~rho~"="~.(sprintf("%1.2f",tempCorvals$estimate))~
               ","~italic(P)~"="~.(sprintf("%1.2f",tempCorvals$p.value))~" )"),
      line=-0.1, family = "Times New Roman", cex=0.65)
mtext(LETTERS[7], 3, outer=FALSE, cex=1.25, adj=-0.21,family = "Times New Roman", line=1.55)
# dev.off()

# svglite(filename = "C:/Users/carac/Dropbox/Analysis/Multivariate_MR/Multi_Latex_Documents/PNGs and EPS/MajorEigentensorEigVecs.svg",
#     width = 11, height = 5.8, pointsize = 10,
#     bg = "white", system_fonts = "Times New Roman")
# grDevices::cairo_ps(filename =
#                       "C:/Users/carac/Dropbox/Analysis/Multivariate_MR/Multi_Latex_Documents/PNGs and EPS/MajorEigentensorEigVecs.eps",
#                     family = "Times", bg="transparent",pointsize=10, width = 11, height = 5.8)

#Plot the major 5 eigentensor eigenvectors
par(mfrow=c(3,5),mar=c(4.5,6,3.5,1)); t=1
for (p in c("M", "R")) {
  for (evec in c(11,16,21,31,36)) {
    if (p=="M"){
      plot(NULL, xlim=c(0.5, 6.5),bty="L",yaxt="n",xlab="Generations",
           ylim=c(-0.05, 1), main="",
           family = "Times New Roman",
           ylab="Among-line varaince \u00B1 90% CI")
      axis(side=2, at=seq(0,1,0.2),family = "Times New Roman", labels = sprintf("%.2f",seq(0,1,0.2)), las=2)
      abline(h=0, lty=2)
      mtext(bquote(italic("e")[.(evec)]),3, adj=0.03, line = 0,family = "Times New Roman")
    } else{
      plot(NULL, xlim=c(0.5, 6.5),bty="L",yaxt="n",xlab="Generations",
           ylim=c(0.5, 2), main="",
           family = "Times New Roman",
           ylab="Residual variance \u00B1 90% CI")
      axis(side=2, at=seq(0.6,2,0.2),family = "Times New Roman", labels = sprintf("%.2f",seq(0.6,2,0.2)), las=2)
      mtext(bquote(italic("e")[.(evec)]),3, adj=0.03, line = 0,family = "Times New Roman")
    } 
    Project_major_EigtenVecs[Parm==p & vec == evec,segments(x,loCI,x, upCI), by=.I]
    Project_major_EigtenVecs[Parm==p & vec == evec,segments(x-0.05,loCI,x+0.05, loCI), by=.I]
    Project_major_EigtenVecs[Parm==p & vec == evec,segments(x-0.05,upCI,x+0.05, upCI), by=.I]
    Project_major_EigtenVecs[Parm==p & vec == evec,
                             points(x,var, pch=21, bg=ifelse(Treat=="1","black", "white")), by=.I]
    mtext(LETTERS[t], 3, outer=FALSE, cex=1.25, adj=-0.45,family = "Times New Roman", line=1.6); t=t+1
  }
}
#look at relationships
for (evec in c(11,16,21,31,36)) {
  Project_major_EigtenVecs_w[vec==evec] %T>% 
    with(plot(R,M, bty="L",  yaxt="n", #xlim=c(0.8, 1.8),ylim=c(0,0.5),
              xlab="Residual variance", ylab="Among-line variance",
              pch=ifelse(Treat==1, 16,21),family = "Times New Roman")) %>% 
    with(text(R,M,Gen,pos=4,cex=0.75, family = "Times New Roman", xpd=TRUE))
  axis(side=2, at=seq(0,0.5, 0.1),family = "Times New Roman", 
       labels = sprintf("%1.1f",seq(0,0.5, 0.1)), las=2)
  mtext(bquote(italic("e")[.(evec)]),3, line = 1.1,family = "Times New Roman")
  met=ifelse(evec %in% c(5),'spearman',"pearson")#16,31
  tempCorvals<-Project_major_EigtenVecs_w[vec==evec, cor.test(R,M, method=met)][c(4,3)]
  if (evec %in% c(5)){#16,31
    mtext(bquote("(Spearman's "~rho~"="~.(sprintf("%1.2f",tempCorvals$estimate))~
                   ","~italic(P)~"="~.(sprintf("%1.2f",tempCorvals$p.value))~" )"),
          line=0, family = "Times New Roman", cex=0.65)
  }else{
    mtext(bquote("(Pearson's "~italic(r)~"="~.(sprintf("%1.2f",tempCorvals$estimate))~
                   ","~italic(P)~"="~.(sprintf("%1.2f",tempCorvals$p.value))~" )"),
          line=0, family = "Times New Roman", cex=0.65) 
  }
  mtext(LETTERS[t], 3, outer=FALSE, cex=1.25, adj=-0.45,family = "Times New Roman", line=1.6); t=t+1
}
# dev.off()
# Project_major_EigtenVecs_w[, shapiro.test(R), by="vec"]

#Table for Supplementary
as.data.frame(apply(matrix(REML_S_eigtenvals[1:(11*6),], ncol=6,  byrow =TRUE), 
      2, function(x) sprintf("%1.3f", x)))->REML_S_eigtenvals_Tab
REML_S_eigtenvals_Tab<-cbind(c(1:11),REML_S_eigtenvals_Tab)
colnames(REML_S_eigtenvals_Tab)<-c("Ek", paste0("e",1:6))
# write.table(REML_S_eigtenvals_Tab, file ="C:/Users/carac/Dropbox/Analysis/Multivariate_MR/Tables/eigentensor_eigenvalues.txt",
#           row.names = FALSE, quote = TRUE, sep="\t")

# # What if we take just run the eigenanalysis on the 10,000 MVN_S?
# MVN_S_eigVals<-t(apply(MVN_S[,,],3, function (x) eigen(x)$values))[,1:11]
# MVN_S_eigVals_CIs<-apply(MVN_S_eigVals, 2, rangeFunc90)
# 
# MVN_S_eigVals_CIs<-rbind(eigen(REML_S)$values[1:11],MVN_S_eigVals_CIs)
# as.data.frame(t(MVN_S_eigVals_CIs))->MVN_S_eigVals_CIs
# colnames(MVN_S_eigVals_CIs)<-c("lam","n","lo", "up")
# setDT(MVN_S_eigVals_CIs)
# par(mfrow=c(1,1))
# plot(NULL, ylim=c(0, 0.06), xlim=c(1,11), ylab="VL",bty="L", xlab="Tensor")
# MVN_S_eigVals_CIs[, points(.I, lam, pch=16), by=.I]
# MVN_S_eigVals_CIs[, segments(.I, lo,.I, up, pch=16), by=.I]
# #Nope, this is wrong
```

```{r Stability in M}
# Eccentricity_CI<-data.frame(NULL)
# for (i in 1:p) {
#   M_vec<-eigen(M_list[[i]])$vector
#   prj_Meigvec<-array(NA, dim=c(MVNsample,6))
#   for (j in 1:6) {
#     prj_Meigvec[,j]<-apply(AsycovM_array[,,,i],3, function (x) t(M_vec[,j]) %*% x %*% t(t(M_vec[,j])))
#   }
#  Eccentricity_CI<- rbind(Eccentricity_CI,
#         c(i,rangeFunc90(prj_Meigvec[,1]/rowSums(prj_Meigvec[,1:6]))))
# }
# colnames(Eccentricity_CI)<-c("p", "n", "loCI", "upCI")
# 
# par(mfrow=c(1,1))
# plot(NULL, ylim=c(0,1), ylab="Eccentricity of M \u00B1 90% CI", 
#      xlim=c(0.5, 6.5),bty="L", xlab="Generations")
# abline(h=1, lty=2)
# for (i in 1:p) {
#   ecc<-(eigen(M_list[[i]])$value[1])/sum(eigen(M_list[[i]])$value[1:6])
#   if (i<=6){
#     points(i-0.1, ecc, pch=16)
#     Eccentricity_CI[Eccentricity_CI$p==i,] %T>% 
#       with(segments(i-0.1, loCI,i-0.1, upCI)) %T>% 
#       with(segments(i-0.15, loCI, i-0.05, loCI)) %>% 
#       with(segments(i-0.15, upCI, i-0.05, upCI))
#   }else{
#     Eccentricity_CI[Eccentricity_CI$p==i,] %T>% 
#       with(segments(i-5.9, loCI,i-5.9, upCI)) %T>% 
#       with(segments(i-5.85, loCI, i-5.95, loCI)) %>% 
#       with(segments(i-5.85, upCI, i-5.95, upCI))
#     points(i-5.9, ecc, pch=21, bg="white")
#   }
# }
```

```{r Vector correlation}
rad2deg <- function(rad) {(rad * 180) / (pi)}
e11<-t(t(REML_S_eigtenvecs[1,])); e16<-t(t(REML_S_eigtenvecs[6,]))*-1
e21<-t(t(REML_S_eigtenvecs[7,]))*-1
e31<-t(t(REML_S_eigtenvecs[13,])); e36<-t(t(REML_S_eigtenvecs[18,]))
h1<- t(t(H_vecs[,1]));h2<- t(t(H_vecs[,2])); h3<- t(t(H_vecs[,3]))*-1

# Dotproduct table
mjrEvecs<-cbind(e11, e16,e21, e31, e36)
DotProd<-rbind(apply(mjrEvecs, 2, function (x) t(h1)%*% x),
               apply(mjrEvecs, 2, function (x) t(h2)%*% x),
               apply(mjrEvecs, 2, function (x) t(h3)%*% x))
# apply(DotProd,2, function (x) round(x,digits=2)) %>% View()

# apply(DotProd,2, function (x) rad2deg(acos(x)))


# mjrVecs<-list(e11, e16, h1,h2,h3)
# mjrVecsNames<-c("e11", "e16", "h1","h2","h3")
# 
# for (i in 1:5) {
#   for (j in 1:5) {
#     cat(mjrVecsNames[i], " and ", mjrVecsNames[j], ":\n",
#       round(rad2deg(acos(t(mjrVecs[[i]]) %*% mjrVecs[[j]])),digits = 0), "\n\n")
#   }
# }
# #S
rad2deg(acos(t(eigen(M_list[[1]])$vector[,1])%*% t(t(eigen(M_list[[2]])$vector[,1]))))
# rad2deg(acos(t(eigen(M_list[[2]])$vector[,1])%*% t(t(eigen(M_list[[3]])$vector[,1]))))
# rad2deg(acos(t(eigen(M_list[[3]])$vector[,1])%*% t(t(eigen(M_list[[4]])$vector[,1]))))
# rad2deg(acos(t(eigen(M_list[[4]])$vector[,1])%*% t(t(eigen(M_list[[5]])$vector[,1]))))
# rad2deg(acos(t(eigen(M_list[[5]])$vector[,1])%*% t(t(eigen(M_list[[6]])$vector[,1]))))
# #L
# rad2deg(acos(t(eigen(M_list[[7]])$vector[,1])%*% t(t(eigen(M_list[[8]])$vector[,1]))))
# rad2deg(acos(t(eigen(M_list[[8]])$vector[,1])%*% t(t(eigen(M_list[[9]])$vector[,1]))))
# rad2deg(acos(t(eigen(M_list[[9]])$vector[,1])%*% t(t(eigen(M_list[[10]])$vector[,1]))))
# rad2deg(acos(t(eigen(M_list[[10]])$vector[,1])%*% t(t(eigen(M_list[[11]])$vector[,1]))))
# rad2deg(acos(t(eigen(M_list[[11]])$vector[,1])%*% t(t(eigen(M_list[[12]])$vector[,1]))))
```


```{r Coordinates analyses}
#############################################################################################################################
#get confidence intervals on the VL Coordinates using the 10,00 REML-MVN estimates of M
# set up Coords array
VLcoordsMVN<-array(NA,dim=c(1,1,m,neigten,MVNsample))

for (i in 1:MVNsample) {
  for (j in 1:neigten) {
    for (k in 1:p){
      frobenius.prod(REML_S_eigTenmat[,,j], AsycovM_array[,,i,k])->VLcoordsMVN[,,k,j,i]
    }
  }
}

m=12
CoordVL_CIs<-as.data.frame(NULL)
for (j in 1:neigten) {
  for (k in 1:p){
    VLcoordsMVN[,,k,j,] %>% rangeFunc90() %>% t %>% as.data.frame %>% 
      magrittr::set_colnames(.,c("n", "lowCI", "upCI")) %>% 
      mutate(pop=k, cord_Eig=paste0("C_E",j)) %>% 
      .[,c(4,5,1:3)] %>% rbind(CoordVL_CIs,.)->CoordVL_CIs
  }
}

#Add reml esitmates
cords_withpercent %>% left_join(CoordVL_CIs, by=c("pop", "cord_Eig")) ->CoordVL_CIs_fin

#max eigentensor number we want to explore
eigtenMAX=6
coordNames<- c(unique(CoordVL_CIs_fin$cord_Eig))[1:eigtenMAX]
CoordVL_CIs_fin %<>% mutate(pop=as.numeric(pop))

ord<-data.frame(pop = 1:12,
                   xnum_coords=c(c(1:6)-0.15, c(7:12)-5.85))
CoordVL_CIs_fin %<>% left_join(ord, by=c("pop"))

#plot CIS
# epsilon=0.9


# setEPS()
# postscript("C:/Users/uqcconra/Dropbox/3_ MR data/Analysis/Univariate_MR/Latex_Documents/VL_Coords.eps", 
#            family = "Times", pointsize=14, width =9, height = 3.5)

# svglite(filename = "C:/Users/carac/Dropbox/Analysis/Multivariate_MR/Multi_Latex_Documents/PNGs and EPS/VL_Coords.svg",
#     width = 7, height = 2.6, pointsize = 12,
#     bg = "white", system_fonts = "Times New Roman")
# grDevices::cairo_ps(filename =
#                       "C:/Users/carac/Dropbox/Analysis/Multivariate_MR/Multi_Latex_Documents/PNGs and EPS/VL_Coords.eps",
#                     family = "Times", bg="transparent",pointsize=12, width = 7, height = 2.6)
# par(mfrow=c(1,3), mar = c(3.5, 5, 1 ,1))
epsilon=0.05
for (i in 1:3) {
  CoordVL_CIs_fin[,c(1:4,9,10,11)] %>% filter(cord_Eig==coordNames[i]) %T>% 
    with(plot(xnum_coords, Coord, bty="L", yaxs="i", main="",family = "Times New Roman", 
              xlab="", xlim=c(0.5,6.5), ylim = c(floor((min(lowCI)*5))/5, ceiling((max(upCI)*5))/5), 
              pch=16, cex=1.2,yaxt="n", xaxt="n", ylab = bquote(paste("Coordinates of "~bold(E)[{.(i)}]~" in"~bold(M)~" \u00B1 90% CI")))) %T>%  
    with(abline(h=0, lty=2)) %T>% 
              #ylab = bquote(paste(italic(C)[bold(M)]^{italic(p)*", "*.(i)})))) %T>% 
    with(., segments(xnum_coords-epsilon, lowCI, xnum_coords+epsilon, lowCI)) %T>% 
    with(., segments(xnum_coords-epsilon, upCI, xnum_coords+epsilon, upCI)) %T>% 
    with(.,segments(xnum_coords, lowCI, xnum_coords, upCI)) %T>% 
    # with(.,axis(1, at=xnum_coords, 
                # labels=as.expression(lapply(pop, function(i)bquote(bold("M")[.(i)]))), cex.axis=.8)) %T>% 
    with(., axis(1, at=seq(1, 6, 1),labels=c(1:6),family = "Times New Roman",
                 tick = TRUE)) %T>%
    with(., axis(1, at=3.5, "Generations",family = "Times New Roman",
                 line=1, tick=FALSE, cex=1.2)) %>% 
    filter(pop %in% 7:12) %>% 
    with(., points(xnum_coords,Coord, pch=21, bg="white"))
  if (i==1){Pseq=seq(-0.4, 0.6, 0.2)}else if(i==2){Pseq=seq(-0.2, 0.6, 0.2)}
  else if(i==4){Pseq=seq(-0.8, 0, 0.2)}  else { Pseq=seq(-0.4, 0.4, 0.2)}
   axis(side=2, at= Pseq, family = "Times New Roman", labels=sprintf("%1.1f",Pseq), las=2) 
  # legend("topright", pch =c(16, 21), bg="white", legend = c("S", "B"), bty = "n")
}
# dev.off()


```

```{r R_eigentensor}
REML_S_R<-array(NA, dim=c(neigten,neigten))
dimnames(REML_S_R)<-list(paste0("e", 1:neigten), paste0("e", 1:neigten))
REML_Rvarmat<-t(apply(R_array, 3, diag))
REML_Rcovmat<-t(apply(R_array, 3, lowerTriangle))
REML_S_R[1:n, 1:n]<- cov(REML_Rvarmat,REML_Rvarmat) #upper left quarter of S
REML_S_R[(n+1):neigten, (n+1):neigten]<-2*cov(REML_Rcovmat,REML_Rcovmat)#lower right quarter of S
REML_S_R[1:n, (n+1):neigten]<-sqrt(2)*cov(REML_Rvarmat, REML_Rcovmat)
REML_S_R[(n+1):neigten, 1:n]<-sqrt(2)*cov(REML_Rcovmat, REML_Rvarmat)
# Eigenanalyses
REML_S_R_eigvec <- eigen(REML_S_R)$vectors
REML_S_R_eigval <- eigen(REML_S_R)$values

#create MVN-sampling S of R
# Get 10,000 S matrices from the 10,000 REML-MVN M matrix sets 
MVN_S_R<-array(NA, dim=c(neigten,neigten,MVNsample))
dimnames(MVN_S_R)<-list(paste0("e", 1:neigten), paste0("e", 1:neigten))
for (i in 1:MVNsample){
  # 1. Construction of M covariance tensor
  MVN_Rvarmat<-t(apply(AsycovR_array[,,i,], 3, diag))
  MVN_Rcovmat<-t(apply(AsycovR_array[,,i,], 3, lowerTriangle))
  MVN_S_R[1:n, 1:n, i]<- cov(MVN_Rvarmat,MVN_Rvarmat) #upper left quarter of S
  MVN_S_R[(n+1):neigten, (n+1):neigten, i]<-2*cov(MVN_Rcovmat,MVN_Rcovmat)#lower right quarter of S
  MVN_S_R[1:n, (n+1):neigten,i]<-sqrt(2)*cov(MVN_Rvarmat, MVN_Rcovmat)
  MVN_S_R[(n+1):neigten, 1:n, i]<-sqrt(2)*cov(MVN_Rcovmat, MVN_Rvarmat)
}

# Get CI for R eigentensors
MVN_S_R_val <- matrix(NA, MVNsample, neigten)
colnames(MVN_S_R_val) <- paste("E", 1:neigten, sep="")
for (i in 1:MVNsample){
  for(j in 1:neigten){
    MVN_S_R_val[i,j] <- t(REML_S_R_eigvec[,j]) %*% MVN_S_R[,,i] %*% REML_S_R_eigvec[,j]
  }
}

# Plot the 90% CIs for eigtensors -----------------------------------------
#Add REML-MVN CI
Eigten_R_CIs<-as.data.frame(NULL)
for (p in 1:neigten){
  Eigten_R_CIs<-rbind(Eigten_R_CIs, c(p, rangeFunc90(MVN_S_R_val[,p])))
}
colnames(Eigten_R_CIs)<-c("Eigten_Num", "n", "lowCI", "upCI")
Eigten_R_CIs %<>% mutate(Eigten=paste0("E", Eigten_R_CIs$Eigten_Num))
# add eigentensor of reml esitmates of R 
Eigten_R_CIs$val<-REML_S_R_eigval
setDT(Eigten_R_CIs)

# svglite(filename = "C:/Users/carac/Dropbox/Analysis/Multivariate_MR/Multi_Latex_Documents/PNGs and EPS/DaveFig_nMR_residual.svg",
#     width = 4, height = 3, pointsize = 10,
#     bg = "white", system_fonts = "Times New Roman")
par(mfrow=c(1,1), mar=c(2,4,2,0.5))
Eigten_R_CIs[1:11, 
             plot(Eigten_Num, val, ylim=c(0,0.14),xlab="",
                  pch=16, bty="L", yaxt="n",xaxt="n", family = "Times New Roman",
                  ylab="Residual variance \u00B1 90% CI")]
axis(side=2, at=seq(0,0.14, 0.02), las=2, family = "Times New Roman",
     labels = sprintf("%1.2f",seq(0,0.14, 0.02)))
Eigten_R_CIs[1:11,segments(Eigten_Num, lowCI,Eigten_Num,upCI)]
Eigten_R_CIs[1:11,segments(Eigten_Num-0.1, lowCI,Eigten_Num+0.1,lowCI)]
Eigten_R_CIs[1:11,segments(Eigten_Num-0.1, upCI,Eigten_Num+0.1,upCI)]
axis(1, at=c(1:Eigtenlim), family = "Times New Roman",
     labels=as.expression(lapply(1:Eigtenlim, function(i)bquote(bold("E")[.(i)]))))
# dev.off()

#Generate the eigentensors matrices
REML_S_R_eigTenmat<-array(NA, dim=c(n,n,neigten))
dimnames(REML_S_R_eigTenmat)<-list(traitnumber, traitnumber,paste0("E", 1:neigten))
for (i in 1:neigten){
  REML_emat<-matrix(NA,n,n)
  lowerTriangle(REML_emat)<-1/sqrt(2)*REML_S_R_eigvec[(n+1):neigten,i]
  REML_emat<-REML_emat+t(REML_emat)
  diag(REML_emat)<-REML_S_R_eigvec[1:n,i]
  REML_S_R_eigTenmat[,,i]<-REML_emat
}

REML_S_R_eigten<-data.frame(NULL)
for (i in 1:neigten){REML_S_R_eigten<-rbind(REML_S_R_eigten,REML_S_R_eigTenmat[,,i])}

REML_S_R_eigtenvecs<-matrix(nrow=n*neigten, ncol=n)
REML_S_R_eigtenvals<-matrix(nrow=n*neigten, ncol=1)
for (i in 1:neigten){ #Using Emma's Method
  REML_S_R_eigtenvecs[((i-1)*n+1):(i*n),]=t(eigen(REML_S_R_eigten[((i-1)*n+1):(i*n),])$vectors) # eigenvectors in rows!!!!!!
  REML_S_R_eigtenvals[((i-1)*n+1):(i*n),1]=eigen(REML_S_R_eigten[((i-1)*n+1):(i*n),])$values
}


```

```{r The two Eigentensors}
# svglite(filename = "C:/Users/carac/Dropbox/Analysis/Multivariate_MR/Multi_Latex_Documents/PNGs and EPS/DaveFigTwoEigtensors.svg",
#     width = 7, height = 3.7, pointsize = 10,
#     bg = "white", system_fonts = "Times New Roman")

# grDevices::cairo_ps(filename =
#                       "C:/Users/carac/Dropbox/Analysis/Multivariate_MR/Multi_Latex_Documents/PNGs and EPS/DaveFigTwoEigtensors.eps",
#                     family = "Times", bg="transparent",pointsize=10, width = 7, height = 3.7)
par(mfrow=c(1,2), mar=c(2,4,2,0.5))
Eigtenlim<-11
plot(NULL, type="n", xlab="", yaxs="i", ylab="Among-line variance \u00B1 90% CI", 
     main="",  xaxt="n", yaxt="n",family = "Times New Roman",
     xlim=c(1,Eigtenlim), ylim=c(0,0.055), bty="L")
epsilon=0.125
for (i in c(1:Eigtenlim)){
  plot_EigtenCIs %>% filter(Eigten_Num==i)->subEig
  with(subEig, segments(Eigten_Num, lowCI, Eigten_Num, upCI))
  with(subEig, segments(Eigten_Num-epsilon, lowCI, Eigten_Num+epsilon, lowCI))
  with(subEig, segments(Eigten_Num-epsilon, upCI, Eigten_Num+epsilon, upCI))
  #with(subEig, points(Eigten_Num, Median))
  with(subEig, points(Eigten_Num, REML_Var, pch=16, cex=0.75))
}
axis(2, at=seq(0,0.055, 0.01), labels=sprintf("%.2f", seq(0,0.055, 0.01)), las=2,family = "Times New Roman")
axis(1, at=c(1:Eigtenlim), family = "Times New Roman",
     labels=as.expression(lapply(1:Eigtenlim, function(i)bquote(bold("E")[.(i)]))))
mtext("A", 3, outer=FALSE, cex=1.25,family = "Times New Roman", adj=-0.25, line=0.65)


Eigten_R_CIs[1:11, 
             plot(Eigten_Num, val, ylim=c(0,0.12),xlab="",yaxs="i",
                  pch=16, bty="L", yaxt="n",xaxt="n", family = "Times New Roman",
                  cex=0.75, ylab="Residual variance \u00B1 90% CI")]
axis(side=2, at=seq(0,0.12, 0.02), las=2, family = "Times New Roman",
     labels = sprintf("%1.2f",seq(0,0.12, 0.02)))
Eigten_R_CIs[1:11,segments(Eigten_Num, lowCI,Eigten_Num,upCI)]
Eigten_R_CIs[1:11,segments(Eigten_Num-epsilon, lowCI,Eigten_Num+epsilon,lowCI)]
Eigten_R_CIs[1:11,segments(Eigten_Num-epsilon, upCI,Eigten_Num+epsilon,upCI)]
axis(1, at=c(1:Eigtenlim), family = "Times New Roman",
     labels=as.expression(lapply(1:Eigtenlim, function(i)bquote(bold("E")[.(i)]))))
mtext("B", 3, outer=FALSE, cex=1.25,family = "Times New Roman", adj=-0.25, line=0.65)
# dev.off()
```



```{r Coordinates of R eigentensors procected onto M}
#get confidence intervals on the R Coordinates using the 10,00 REML-MVN estimates of M
# set up Coords array
VLthruRcoordsMVN<-array(NA,dim=c(1,1,m,neigten,MVNsample))

for (i in 1:MVNsample) {
  for (j in 1:neigten) {
    for (k in 1:p){
      frobenius.prod(REML_S_eigTenmat[,,j], AsycovR_array[,,i,k])->VLthruRcoordsMVN[,,k,j,i]
    }
  }
}

m=12
CoordVLthruR_CIs<-as.data.frame(NULL)
for (j in 1:neigten) {
  for (k in 1:p){
    VLthruRcoordsMVN[,,k,j,] %>% rangeFunc90() %>% t %>% as.data.frame %>% 
      magrittr::set_colnames(.,c("n", "lowCI", "upCI")) %>% 
      mutate(pop=k, cord_Eig=paste0("C_E",j)) %>% 
      .[,c(4,5,1:3)] %>% rbind(CoordVLthruR_CIs,.)->CoordVLthruR_CIs
  }
}
setDT(CoordVLthruR_CIs)
CoordVLthruR_CIs[,EigNum:=as.numeric(gsub("C_E", "",cord_Eig))]
CoordVLthruR_CIs[,p:=as.numeric(pop)]
CoordVLthruR_CIs[,Coord:=0]
for (i in 1:neigten) {
  for (j in 1:p) {
    CoordVLthruR_CIs[EigNum==i & p==j]$Coord<-frobenius.prod(REML_S_eigTenmat[,,i], R_list[[j]])
  }
}


#max eigentensor number we want to explore
eigtenMAX=6
coordNames<- c(unique(CoordVLthruR_CIs$cord_Eig))[1:eigtenMAX]
ord<-data.frame(pop = 1:12,
                   xnum_coords=c(c(1:6)-0.15, c(7:12)-5.85))
CoordVLthruR_CIs %<>% left_join(ord, by=c("pop"))

#plot CIS
# epsilon=0.9


# setEPS()
# postscript("C:/Users/uqcconra/Dropbox/3_ MR data/Analysis/Univariate_MR/Latex_Documents/VL_Coords.eps", 
#            family = "Times", pointsize=14, width =9, height = 3.5)

# svglite(filename = "C:/Users/carac/Dropbox/Analysis/Multivariate_MR/Multi_Latex_Documents/PNGs and EPS/VLthruR_Coords.svg",
#     width = 7, height = 2.8, pointsize = 12,
#     bg = "white", system_fonts = "Times New Roman")
# grDevices::cairo_ps(filename =
#                       "C:/Users/carac/Dropbox/Analysis/Multivariate_MR/Multi_Latex_Documents/PNGs and EPS/VLthruR_Coords.eps",
#                     family = "Times", bg="transparent",pointsize=12, width = 7, height = 2.6)
par(mfrow=c(1,3), mar = c(3.5, 5, 1 ,1))
par(mfrow=c(1,3), mar = c(3.5, 5, 1 ,1))
epsilon=0.05
for (i in 1:3) {
  CoordVLthruR_CIs[cord_Eig==coordNames[i]] %T>% 
    with(plot(xnum_coords, Coord, bty="L", yaxs="i", main="",family = "Times New Roman", 
              xlab="", xlim=c(0.5,6.5), ylim = c(floor((min(lowCI)*5))/5, ceiling((max(upCI)*5))/5), 
              pch=16, cex=1.2,yaxt="n", xaxt="n", ylab = bquote(paste("Coordinates of "~bold(E)[{.(i)}]~"in"~bold(R)~" \u00B1 90% CI")))) %T>%  
              #ylab = bquote(paste(italic(C)[bold(M)]^{italic(p)*", "*.(i)})))) %T>% 
    with(., segments(xnum_coords-epsilon, lowCI, xnum_coords+epsilon, lowCI)) %T>% 
    with(., segments(xnum_coords-epsilon, upCI, xnum_coords+epsilon, upCI)) %T>% 
    with(.,segments(xnum_coords, lowCI, xnum_coords, upCI)) %T>% 
    # with(.,axis(1, at=xnum_coords, 
                # labels=as.expression(lapply(pop, function(i)bquote(bold("M")[.(i)]))), cex.axis=.8)) %T>% 
    with(., axis(1, at=seq(1, 6, 1),labels=c(1:6),family = "Times New Roman",
                 tick = TRUE)) %T>%
    with(., axis(1, at=3.5, "Generations",family = "Times New Roman",
                 line=1, tick=FALSE, cex=1.2)) %>% 
    filter(pop %in% 7:12) %>% 
    with(., points(xnum_coords,Coord, pch=21, bg="white"))
  Pseq=seq(-1.8, 1.6, 0.2)
  axis(side=2, at= Pseq, family = "Times New Roman", labels=sprintf("%1.1f",Pseq), las=2) 
  # legend("topright", pch =c(16, 21), bg="white", legend = c("S", "B"), bty = "n")
}
# dev.off()


```


```{r Mitteroecker Distance}
# pr.coord (principle co-ordinates analysis)
# doesn't work for R 
# doesn't work for Eigvectors as one M is not positive definite
library(vcvComp)
citation("vcvComp")

Hspace_arrayM<-array(NA, dim=c(3,3,12))
for(i in 1:12){
  Hspace_arrayM[,,i]<-t(H_vecs) %*% M_array[,,i] %*% H_vecs
}

# check matrices are pos.def
apply(Hspace_arrayM[,,],3,is.positive.definite)


# Now, to calculate the D (distance matrix) and get the principle 
disMAtH_M<-matrix(NA, nrow = 12, ncol = 12)
for (i in 1:12){
  k=i
  while (k<=12) {
    Hspace_arrayM[,,i] %*% solve(Hspace_arrayM[,,k])->tmp
    sqrt(sum(log(eigen(tmp)$values)^2))->disMAtH_M[i,k]
    k=k+1
  }
}
lowerTriangle(disMAtH_M) = upperTriangle(disMAtH_M, byrow=TRUE)

# pr.coord (principle co-ordinates analysis)
vcvComp::pr.coord(disMAtH_M)->prcoordSupH_M


# Scree plot of principle coordinates
plot_prcoordSupH<-as.data.frame(cbind(1:11,prcoordSupH_M$values[1:11]))
colnames(plot_prcoordSupH)<-c("Cord","Eigval")
setDT(plot_prcoordSupH)
plot_prcoordSupH[,Eigval:=as.numeric(Eigval)]
plot_prcoordSupH[,Prop:=round(Eigval/sum(Eigval)*100, 1)]

as.data.frame(cbind(1:12, prcoordSupH_M$PCoords[,1:3]))->PlotPcords
colnames(PlotPcords)[1]<-c("p")
PlotPcords$pop=unlist(names(M_list))
setDT(PlotPcords)
PlotPcords[, c("Treat","Gen") := data.table(str_split_fixed(pop,"_",2))]
PlotPcords[, c("Treat","Gen") := PlotPcords[, lapply(.SD, as.numeric), .SDcols =c("Treat","Gen")]]

#plot
# svglite(filename = "C:/Users/carac/Dropbox/Analysis/Multivariate_MR/Multi_Latex_Documents/PNGs and EPS/M_in_H_MittDist.svg",
#     width = 7.5, height = 5.4, pointsize = 11,
#     bg = "white", system_fonts = "Times New Roman")
# grDevices::cairo_pdf(filename =
#                       "C:/Users/carac/Dropbox/Analysis/Multivariate_MR/Multi_Latex_Documents/PNGs and EPS/M_in_H_MittDist.pdf",fallback_resolution = 2400,
#                     family = "Times", bg="transparent",pointsize=11, width = 7.5, height = 5.4)

#including a 3d plot
# grDevices::cairo_pdf(filename =
#                       "C:/Users/carac/Dropbox/Analysis/Multivariate_MR/Multi_Latex_Documents/PNGs and EPS/M_in_H_MittDist_3D.pdf",fallback_resolution = 2400,
#                     family = "Times", bg="transparent",pointsize=11, width = 7.25, height = 8.4)
# svglite(filename = "C:/Users/carac/Dropbox/Analysis/Multivariate_MR/Multi_Latex_Documents/PNGs and EPS/M_in_H_MittDist_3D.svg",
#     width = 7.25, height = 8.4, pointsize = 11,
#     bg = "white", system_fonts = "Times New Roman")
m<-matrix(c(rep(1:2, each = 3),
  rep(1:2, each = 3),
  rep(3:4, each = 3),
  rep(3:4, each = 3),
  rep(5, 18)), 7,6, byrow = TRUE)
layout(m)
par(mar = c(5,5,2,1))# mfrow=c(2,2),
plot_prcoordSupH[,plot(Cord, Prop, bty="L",ylim=c(0,40),
                       ylab="", xlab="Principal coordinate",
                       xlim=c(0,12),xaxs="i",family = "Times New Roman",
                       yaxs="i",pch=16, xaxt="n", yaxt="n")]
axis(1, at=c(seq(1,11,1)), labels=seq(1,11,1),family = "Times New Roman", cex.axis=1)
axis(2, at=c(seq(0,40, 5)), las=2,cex=0.6, family = "Times New Roman",
     labels=paste0(format(round(seq(0,40, 5),1), nsmall = 0),"%"))
axis(2, at=20, labels = "Proportion of varaince",family = "Times New Roman", line=2.1, tick=F)
# plot_prcoordSupH[,text(Cord, Prop, Prop, pos=3, cex=0.8)]
mtext("A", 3,family = "Times New Roman", outer=FALSE, cex=1.35,adj=-0.22, line=.65)
# 47% of the information captured in the first 2 principle coordinates, and 58% in 3
#  2d with lines
par(xpd=TRUE, family = "Times New Roman")
for (i in 1:3) {
  if (i==1){
    Cols<-c("PCo1","PCo2")
    xlimy=c(-1,2);ylimy=c(-1,1)
    lab<-c("First principal coordinate", "Second principal coordinate")
    labpos_1=c(2,2,1,1,2,4)
    labpos_3=c(4,2,1,1,4,2)
  }else if (i==2){
    Cols<-c("PCo1","PCo3")
    xlimy=c(-1,2);ylimy=c(-1.5,1)
    lab<-c("First principal coordinate", "Third principal coordinate")
    labpos_1=c(4,2,4,4,3,4)
    labpos_3=c(4,3,1,3,1,1)
  }else {
    Cols<-c("PCo2","PCo3")
    xlimy=c(-1,1.5);ylimy=c(-1.5,1)
    lab<-c("Second principal coordinate", "Third principal coordinate")
    labpos_1=c(1,2,1,2,3,2)
    labpos_3=c(4,2,2,1,3,4)
  }
  plot(NULL, ylim=ylimy, yaxt="n", ylab=lab[2],family = "Times New Roman",
       xlab=lab[1], bty="L", xlim=xlimy)
  axis(side=2, at=seq(ylimy[1],ylimy[2],0.5), las=2, family = "Times New Roman",
       labels=sprintf("%1.1f",seq(ylimy[1],ylimy[2],0.5)))
  abline(h=0, lty=2, xpd=F);abline(v=0, lty=2, xpd=F)
  PlotPcords[Treat==1, lines(.SD, lwd=2, 
                             col=rgb(0.2,0.2,0.2,alpha=0.75)),.SDcols = Cols]
  PlotPcords[Treat==1, points(.SD, pch=16, cex=1.15),.SDcols = Cols]
  PlotPcords[Treat==3, lines(.SD, lwd=2, col=rgb(0.1,0.2,1,alpha=0.75)),.SDcols = Cols]
  PlotPcords[Treat==3, points(.SD, pch=16, cex=1.15, col="blue"),.SDcols = Cols]
  PlotPcords[Treat==1,text(.SD, labels=Gen, family = "Times New Roman",pos=labpos_1),.SDcols = Cols]
  PlotPcords[Treat==3,text(.SD, labels=Gen, family = "Times New Roman",pos=labpos_3, col="blue"),.SDcols = Cols]
  mtext(LETTERS[i+1], 3,family = "Times New Roman", outer=FALSE, cex=1.35,adj=-0.22, line=.65)
}
legend("bottomright",pch=16, col=c(1,"blue"),lty=1,
       c("Small population", "Large Population"),bty="n")
# dev.off()

#Offset labels
par(mar = c(2,5,2,1))
PC_off<-PlotPcords[,.(p, pop, Treat, Gen)]
PC_off$x=c(0,-0.15,0,-0.15,0,0.06, 0.08,-0.12,0,0,0.07,0)
PC_off$z=c(-0.09,0,-0.11,0,0.09,0, 0,0,-0.09,-0.07,0,0.08)

library("scatterplot3d")
PlotPcords_3D<-merge(PlotPcords, PC_off)
PlotPcords_3D %>% as.data.frame %>% 
  with(scatterplot3d(x=PCo1,y=PCo2,z=PCo3,
                            type = "h",pch=16,
                            xlab="First principal coordinate", 
                            ylab="", color=ifelse(Treat==1, "black", "blue"),
                            zlab="Third principal coordinate"))
text(6.5, -2, "Second principal coordinate", srt=50)

library("plot3D")
PlotPcords_3D<-merge(PlotPcords, PC_off)
PlotPcords_3D[Treat==1,scatter3D(x=PCo1,y=PCo2,z=PCo3,pch=16,nticks =6, label=TRUE,
                     col="black", type="h",colkey = FALSE, ticktype="detailed", 
                     axes=TRUE,theta=30, phi=25,bty = "b2",
                     xlim=c(-1, 1.8), ylim=c(-1, 1),zlim=c(-1.1, 0.6),
                     xlab="First principal coordinate", 
                     ylab="Second principal coordinate", 
                     zlab="Third principal coordinate") ]
PlotPcords_3D[Treat==1,scatter3D(x=PCo1,y=PCo2,z=PCo3,pch=16, label=TRUE,
                     col="black",add = TRUE, type="b", lwd=2)]
PlotPcords_3D[Treat==3,scatter3D(x=PCo1,y=PCo2,z=PCo3,pch=16, label=TRUE,
                     col="blue", add = TRUE, type="h")]
PlotPcords_3D[Treat==3,scatter3D(x=PCo1,y=PCo2,z=PCo3,pch=16, label=TRUE,
                     col="blue",add = TRUE, type="b", lwd=2)]
PlotPcords_3D[,text3D(x=PCo1+x,y=PCo2,z=PCo3+z, colkey = FALSE, add = TRUE, 
        labels = Gen, colvar=Treat, col = c("black", "blue"))]
legend(0.5,0,pch=16, col=c(1,"blue"),lty=1,
       c("Small population", "Large Population"),bty="n")
mtext(LETTERS[5], side=3, cex=1.35, adj=0.25,family = "Times New Roman",line=-0.35)
# dev.off()

# ## Original 
# library("plot3D")
# PlotPcords_3D<-merge(PlotPcords, PC_off)
# PlotPcords_3D[Treat==1,scatter3D(x=PCo1,y=PCo2,z=PCo3,pch=16,nticks =6, label=TRUE,
#                      col="black", type="b",colkey = FALSE, ticktype="detailed", 
#                      axes=TRUE,theta=30, phi=25,bty = "u",
#                      xlim=c(-1, 1.5), ylim=c(-1, 1),zlim=c(-1, 0.6),
#                      xlab="First principal coordinate", 
#                      ylab="Second principal coordinate", 
#                      zlab="Third principal coordinate") ]
# PlotPcords_3D[Treat==3,scatter3D(x=PCo1,y=PCo2,z=PCo3,pch=16, label=TRUE,
#                      col="blue", add = TRUE, type="b")]
# PlotPcords_3D[,text3D(x=PCo1+x,y=PCo2,z=PCo3+z, colkey = FALSE, add = TRUE, 
#         labels = Gen, colvar=Treat, col = c("black", "blue"))]
# legend(0.5,0,pch=16, col=c(1,"blue"),lty=1,
#        c("Small population", "Large Population"),bty="n")
# mtext(LETTERS[5], side=3, cex=1.35, adj=0.25,family = "Times New Roman",line=-0.35)
```





# Multivariate Analyses of Mutation

#packages
library(dplyr, warn.conflicts = FALSE); options(dplyr.summarise.inform = FALSE); library(data.table)
library(stringi); library(stringr); library(magrittr); library(evolqg); library(Matrix)
library(foreach); library(matrixStats); library(MASS); library(parallel); library(matrixcalc)
library(gdata); library(bigstatsr); library(abind); library(psych); library(stats)

# # font for eps figure outputs
# library(sysfonts); library(showtextdb); library(showtext)
# ## Load Times New Roman fonts on Windows
# font_add("times", regular = "times.ttf", bold = "timesbd.ttf", italic = "timesi.ttf", bolditalic ="timesbi.ttf")

# library(showtext)
# ## add the Arial font
# font_add("Arial", regular = "arial.ttf",
#          bold = "arialbd.ttf", italic = "ariali.ttf", bolditalic = "arialbi.ttf")

# font for png
library(extrafont)
loadfonts(device = "win")

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
nullnumber <- 1000 # number of randomised datasets
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

theta_array <- array(NA, dim = c((n*(n+1)/2), 1, p))
for (i in c(1:p)){theta_array[,,i] <- theta_list[[i]]}
dimnames(theta_array)[[3]] <- names(theta_list)

# Create array for 10,000 REML-MVN data
AsycovM_array <- array(NA, dim = c(n, n, MVNsample, p)) 

# The loop for REML-MVN 
set.seed(42)
for (i in c(1:p)){ 
  reps <- mvrnorm(n = MVNsample, theta_array[,,i], V_array[,,i])
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
  m <- matrix(NA, nrow = n, ncol = n)
  lowerTriangle(m, diag = TRUE, byrow = TRUE) <- x
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

null_M_array <- foreach (i = 1:length(null_Mfilenames), 
                         .packages = c('data.table', 'foreach','gdata'), 
                         .combine='acomb', .multicombine=TRUE) %dopar% {
                           rdm_null_array <- rando_M_import(i)
                           rdm_null_array
                         }

parallel::stopCluster(cluster)

dim(null_M_array)

# Omit models that haven't converged
for (i in 1:nrow(unconverged_null)) {
  null_M_array[,,unconverged_null$RandRep[i],
               paste(unconverged_null[i, 2:3], collapse = "_")] <- matrix(NA, ncol=6, nrow = 6)
}


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


# Section 2. M matrices eigenanalyses and (co)variance investigation ------------------------------------------------------
# Variability in statistical support & sign
# Output the M matrices and eigenanalyses tables with 90% CI for supplementary
# Create a figure contrasting M matrices across the 12 populations against univariate estimates

# Generate the coordinates of the 6 x 6 matrix
comb <- data.frame(NULL)
l = 0; e = c(1:n)
while (l < 6) {
  comb <- rbind(comb, cbind(rep((l+1), (n-l)), e[(l+1):n]))
  l = l + 1
} 

# Calculated the CI for each element in M using the 10,000 REML-MVN AsycovM_array
# Want 90% for variances, and 95% covariances 
M_tab <- data.frame(NULL)
for (mat in 1:p) {
  for (j in 1:nrow(comb)) {
    vl<-sprintf("%.3f",M_list[[mat]][comb[j,1],comb[j,2]])
    if (comb[j,1] != comb[j,2]){
      vl_ci <- paste(sprintf("%.3f",rangeFunc95(AsycovM_array[comb[j,1], comb[j,2],,mat])[2:3]), collapse="; ")
    } else {
      vl_ci <- paste(sprintf("%.3f",rangeFunc90(AsycovM_array[comb[j,1], comb[j,2],,mat])[2:3]), collapse="; ")
    }
    M_tab <- rbind(M_tab, cbind(names(M_list)[[mat]], comb[j,1], comb[j,2], vl,vl_ci))
  }
}
setDT(M_tab)
M_tab_lng <- M_tab # for stacked barchart
M_tab[, c("Treat","Gen") := data.table(str_split_fixed(V1,"_",2))]
M_tab <- melt(M_tab[,-1], id.var = c("Treat", "Gen", "V2", "V3"), measure.vars = c("vl","vl_ci") )
M_tab <- dcast(M_tab, Gen + variable + V3 ~ Treat + V2, value.var = c("value"))
setorderv(M_tab, c("Gen", "V3"), c(1,1))
M_tab[is.na(M_tab)] <- ""

# write.table(M_tab, file =paste(".",outdir_tab,"M_tab_CI.txt", sep="/"),
#           row.names = FALSE, quote = TRUE, sep="\t")


# Create a stacked bargraph for the covariances
M_tab_lng[, "ci_lo":= unlist(strsplit(vl_ci,"; "))[1], by=.I]
M_tab_lng[, "ci_up":= unlist(strsplit(vl_ci,"; "))[2], by=.I]
colnames(M_tab_lng)[c(1,3)] <- c("pop", "trait_num")
M_tab_lng <- cbind(M_tab_lng[,c(1,5)],apply(M_tab_lng[,c(2:4,6:9)],2, as.numeric))
# Create a trait name dataframe to add back the wing ILDs
Trait_name = data.frame(trait_num = 1:6, 
                        trait = c("CS","1.2","1.5", "2.5","2.8", "3.7"))
M_tab_lng <- merge(M_tab_lng, Trait_name, by= "trait_num")
colnames(M_tab_lng)[c(1, 4,10)] <- c("trait2_num","trait_num","trait2")
M_tab_lng <- merge(M_tab_lng, Trait_name, by= "trait_num")
M_tab_lng[, Sig:= fcase(ci_lo < 0 & ci_up < 0 & vl < 0, 1,
                        ci_lo > 0 & ci_up > 0 & vl > 0, 1,
                        default = 0)]

# write.csv(M_tab_lng, file =paste(".",outdir_tab,"M_tab_lng_sig.csv", sep="/"),
#           row.names = FALSE)

# Now to count the significant covariances and variances
M_cov_count <- M_tab_lng
M_cov_count[, Sign:= fcase(vl >= 0, "pos",
                           default = "neg") ]
M_cov_count <- M_cov_count[, .(N = .N), by=c("trait_num","trait2_num", "Sign", "Sig")]
M_cov_count[,.(N = sum(N)), by=c("trait_num","trait2_num")] # Check that N sums to 12 
M_cov_count[, tn := paste(trait_num,trait2_num, sep="_")]

# Tidy trait names
Trait_name_lng <- unique(M_tab_lng[,.(trait_num,trait2_num,trait, trait2)])
Trait_name_lng[, cov:=fcase(trait_num == trait2_num, "var", default = "cov")]
Trait_name_lng[, tn := paste(trait_num,trait2_num, sep="_")]
Trait_name_lng[, Name := paste(trait,trait2, sep=" - ")]
Trait_name_lng[cov == "var", Name := trait]
setorderv(Trait_name_lng, c("cov","trait_num", "trait2_num"), c(-1,-1,-1))
Trait_name_lng[, x := .I]

# Merge data frames and plot correlation count
M_cov_count_plot <- merge(M_cov_count[,.(Sign,Sig,N,tn)],Trait_name_lng[,.(tn,Name,x, cov)],all.x = TRUE, by="tn")
M_cov_count_plot$Name = factor(M_cov_count_plot$Name, 
                               levels = c("3.7","2.8 - 3.7","2.8","2.5 - 3.7","2.5 - 2.8",
                                          "2.5","1.5 - 3.7","1.5 - 2.8","1.5 - 2.5","1.5",
                                          "1.2 - 3.7", "1.2 - 2.8", "1.2 - 2.5", "1.2 - 1.5", 
                                          "1.2","CS - 3.7","CS - 2.8","CS - 2.5","CS - 1.5",
                                          "CS - 1.2", "CS"))
m_cov_list <- list(); t=1
for (i in c("cov","var")) {
  M_cov_dat <- dcast(M_cov_count_plot[cov==i,.(Sig, Sign, Name, N)], 
                     Sig+Sign~ Name, value.var = "N")
  M_cov_dat[is.na(M_cov_dat)] <- 0
  
  # Get positive values first
  M_cov_dat[, Sig := fcase(Sig== 0, "NS", default = "Sig")]
  M_cov_dat_pos <- as.matrix(M_cov_dat[Sign=="pos",-c(1,2)])
  rownames(M_cov_dat_pos) <- M_cov_dat[Sign=="pos",]$Sig
  m_cov_list[[t]] <- M_cov_dat_pos; t = t + 1
  
  if(i=="cov"){
  # then negative values (varainces will be empty)
  M_cov_dat_neg <- as.matrix(M_cov_dat[Sign=="neg",-c(1,2)])
  rownames(M_cov_dat_neg) <- M_cov_dat[Sign=="neg"]$Sig
  M_cov_dat_neg <- M_cov_dat_neg * -1
  m_cov_list[[t]] <- M_cov_dat_neg; t = t + 1
  } else {
    break
  }
}

# postscript(paste(".",outdir_fig,"m_cov_count.eps", sep="/"),
#            horizontal = FALSE, onefile = FALSE, paper = "special",
#            pointsize = 16.95, width = 6.511811, height = 7.161417)
# # To convert, use png mm width converted into inches / 10; use "Times" for Times New Roman
# 
# png(paste(".",outdir_fig,"m_cov_count.png", sep="/"),  res = 300,
#     bg = "white", pointsize = 14, width = 1654, height = 1819)

# Variance on the top and covariance at the bottom
layout(mat = matrix(c(0,2,1,2,1,2,0,2), 
                    nrow = 2, ncol = 4),
       heights = c(0.333, 0.666))
# layout.show(2)


# Variances first
par(mar=c(4,6,2,1))
barplot(m_cov_list[[3]], 
        col = c("grey", 1),xaxt="n", 
        border = NA,
        horiz = T, xlim=c(0,12),
        xlab = "Count",
        las=2, space = 0.15)
axis(side=2, at= 3.45, labels = "Variance", 
     line = 2, tick = FALSE)
axis(side=1, at = seq(0,12,2),  
     labels = seq(0,12,2))
legend("top",inset = c(0, -0.34), fill= c("grey", 1), bty="n", 
       border = c("grey", 1),
       pt.cex=5,horiz = TRUE,xpd = TRUE, cex=1.05,
       legend = c("Not Significant", "Significant"))
mtext("(a)", 3,outer=FALSE, cex = 1, adj = -0.45, line=1)

par(mar=c(4,6,2,1))
barplot(m_cov_list[[1]], 
        col = c("grey", 1),xaxt="n", 
        border = NA,
        horiz = T, xlim=c(-12,12),
        xlab = "Count",
        las=2, space = 0.15)
axis(side=1, at = seq(-12, 12, 2),  
     labels = seq(-12, 12, 2))
axis(side=2, at= 8.625, labels = "Covariance", 
     line = 3.75, tick = FALSE)
barplot(m_cov_list[[2]], 
        col = c("grey", 1), horiz = T,
        border = NA,
        xaxt="n",
        las=2, xlim=c(-12,12),
        space = 0.15, add=T)
abline(v=0, col="black", lwd = 2.65, lty=2)
mtext("(b)", 3,cex = 1, adj = -0.21)
# dev.off()



# Get CIs for M Eigenvectors
n_cores <- detectCores()-1
cluster <- makeCluster(n_cores)
doParallel::registerDoParallel(cluster)
acomb <- function(...) abind(..., along=4)


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

# Generate a table of all of the Eigenvalues of M
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

# Generate a table of all of the Eigenvectors of M
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

# Generate the figure comparing eigenvalues to individual variances
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

# postscript(paste(".",outdir_fig,"m_eigenvalues.eps", sep="/"),
#            horizontal = FALSE, onefile = FALSE, paper = "special",
#            pointsize = 14.95, width = 4.185, height = 4.7874)

# png(paste(".",outdir_fig,"m_eigenvalues.png", sep="/"),  res=300,
#     bg = "white", pointsize = 10, width = 1063, height = 1036)

par(mfrow=c(1,1), mar=c(4, 4.5,2, 1))
plot(NULL, xlim=c(0.5,6.5), xlab="",yaxt="n", bty="L",
     ylab="",xaxt="n",yaxs="i", ylim=c(0,100))
M_VarNeigenVec[,boxplot(Eig~VecNum,col="white", outcex = 0.75, at=seq(0.85,5.85, 1), whisklty = 1,
                        yaxt="n",pch=4,xaxt="n",frame=F,boxwex=0.2, add = TRUE)]
M_VarNeigenVec[,boxplot(Var~VecNum,col="darkgrey", outcex = 0.75, at=seq(1.15,6.15, 1), whisklty = 1,
                        yaxt="n",pch=4,xaxt="n",frame=F,boxwex=0.2, add = TRUE)]
axis(side=1,at=3.5, "Ordered eigenvalue or trait", line=2, tick = FALSE, cex.axis = 1.2)
axis(side=2,at=50, "Proportion of among-line variance", 
     cex.axis = 1.2, line=2.5, tick = FALSE)
axis(side=2, at=seq(0,100, 20),paste0(seq(0,100, 20),"%"), las=2)
axis(side=1, at=1:6, as.expression(lapply(paste0(1:6,c("st","nd","rd","th","th","th")), function(i)bquote(italic(.(i))))))
# dev.off()




# Section 3.A Krzanowksi's Common Subspaces, H ------------------
# Part A. Implement Krzanowksi's method, and generate eigenanalysis and null distributions 

# We first determined the number of eigenvalues capturing at least 90% variance
M_90var_eigs <- foreach(i = 1:p, .combine = 'rbind') %do%
  {cbind(i, 1:6, 
         eigen(M_list[[i]])$values / tr(M_list[[i]]))}
colnames(M_90var_eigs) <- c("pop", "eigval_num", "perc_variance")
M_90var_eigs <- as.data.table(M_90var_eigs)
M_90var_eigs[, cum_sum := cumsum(perc_variance), by="pop"]
M_90var_eigs <- M_90var_eigs[cum_sum >= .9, .SD[1], pop]
# save the number of eigenvalues (j in equation 2) as a vector
m_90_neigvec <- M_90var_eigs$eigval_num
# make sure that the order is the same as the 12 pops in M_list

# H matrix determined with subset M with eigenvectors capturing at least 90% variance
# takes a list of matrices and vector of # eigenvectors (j estimated above)
# Returns H
# adapted from Melo's evolq R package
KrzSubspace_90var <- function(mat_list, num_eigvec){
  LL_T <- mapply(function(x, y){
    L = eigen(x)$vectors[, 1:y]
    L %*% t(L)
  }, mat_list, num_eigvec, 
  SIMPLIFY = FALSE)
  H = Reduce("+", LL_T)
  H
}

# using the first 3 largest eigenvectors from M to determine in H
# then saving the eigenvectors of H
round(KrzSubspace_90var(M_list, m_90_neigvec), digits = 3) # Inspect the H matrix
H_vecs <-  eigen(KrzSubspace_90var(M_list, m_90_neigvec))$vector[,1:3]
# Flipping vectors, so that the largest loading for h1 and h3 is positive
H_vecs <- H_vecs %*% diag(c(-1,1,-1))
H_val <-  eigen(KrzSubspace_90var(M_list, m_90_neigvec))$values[1:3]

n_cores <- detectCores()-1
cluster <- makeCluster(n_cores)
doParallel::registerDoParallel(cluster)
acomb <- function(...) abind(..., along=3)

remlmvn_H_array <- foreach (i = 1:MVNsample, 
                         .packages = c('foreach'), 
                         .combine='acomb', .multicombine=TRUE) %dopar% {
                           temp_rml_list <- asplit(AsycovM_array[,,i,], 3)
                           KrzSubspace_90var(temp_rml_list, m_90_neigvec)
                         }

null_H_array <- foreach(i = 1:nullnumber, 
                        .packages = c('evolqg'), 
                        .combine='acomb', .multicombine=TRUE) %dopar% {
                          temp_null_list <- asplit(null_M_array[,,i,], 3)
                          temp_null_list <- Filter(function(a) any(!is.na(a)), temp_null_list)
                          KrzSubspace_90var(temp_null_list, rep(3,p))
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
H_eigval_CI <- cbind(1:3, H_val, 
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
setorderv(H_tab, cols = c("Est", "variable"),order = c(-1L, 1L))

# write.csv(H_tab, file =paste(".",outdir_tab,"EigH_tab.csv", sep="/"),
# row.names = FALSE, quote = TRUE)


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


H_thru_M <- cbind(rep(1:12, 6), rep(1:3,each=6), 
                  rep(c(seq(1.25,11.25,2), seq(1.75,11.75,2)),3))
H_thru_M <- as.data.table(H_thru_M)
colnames(H_thru_M)<-c("pop_num", "Hvec","x")
H_thru_M[, pop:=c(unlist(names(M_list)[p])), by=.I]
H_thru_M[,c("var","Lo", "Hi"):=0]

for(j in c(1:3)){#number of eigenvecotrs of H
  for (i in c(1:p)) {
    H_thru_M[pop_num ==i & Hvec==j, var:= Var_HeigVecs[,i,j]]
    tempCI=rangeFunc90(AsycovHeigVecs_array[,,,i,j])
    H_thru_M[ pop_num ==i & Hvec==j, `:=`(Lo=tempCI[2] , Hi=tempCI[3])]
  }
};rm(tempCI)
H_thru_M[,pchy:=fcase(pop_num %in% 1:6, 16,
                        default= 21)]

# postscript(paste(".",outdir_fig,"HvecsThruM.eps", sep="/"),
#            family = "Times", pointsize=12, width =10, height = 4.35)
# 
# png(paste(".",outdir_fig,"HvecsThruM.png", sep="/"), res=300,
#     bg = "white", pointsize = 14, width = 797.25, height = 2162.3)

par(mfrow = c(3,1), mar = c(3,3,2,1), oma = c(2,2,0,0))
for(j in c(1:3)){
  yaxlim = if(j==1){c(0, 0.8)} else if (j ==2) {c(-0.1, 0.7)} else {c(-0.1, 0.3)}
  H_thru_M[Hvec == j,
           plot(x, var, ylim = yaxlim, las=2, xlim = c(0.5,12.5), xaxt = "n",
                ylab = "", xlab = "",
                main = bquote(italic("h")[.(j)]))]
  abline(h=0, lty=2)
  axis(side=1, at=seq(1.5,11.5,2), labels = 1:6)
  H_thru_M[Hvec == j, segments(x,Lo,x, Hi)]
  H_thru_M[Hvec == j, segments((x-0.15), Lo,(x+0.15), Lo)]
  H_thru_M[Hvec == j, segments((x-0.15), Hi,(x+0.15), Hi)]
  H_thru_M[Hvec == j, points(x,var, pch=pchy, bg=ifelse(pchy==21, "white", NA))]
  if (j==1){
    legend("topright", pch =c(16, 21), col = c(1,1),
           bg =c(1,0),
           legend= c("Small", "Large"), bty="n")
  }
}
mtext("Among-line variance \u00B1 90% CI", side = 2.25, cex=0.85, outer = T)
mtext("Generations", side = 1, at = 6.5, line = 3.25,  cex = 0.85)
# dev.off()

# Section 3.B Krzanowksi's Common Subspaces, H  h scores--------------------------
PC_scores <- data %*% PCA$vectors ## calculate PC scores for all objects and PCs
colnames(PC_scores)=c("PC1", "PC2")
head(PC_scores)

# Import data
wing_dat <- fread("mr_wings_6traits.csv")
wing_dat <- wing_dat[,c(1:7, 10, 13)] # using the standardised score
wing_data_wide <- dcast(wing_dat, Animal+Gen+Treatment+Treat+Line+Vial+Ind ~ Index,  
                        value.var ="stdScore")
h_scores <- as.matrix(wing_data_wide[,c(8:13)]) %*% H_vecs
h_scores <- as.data.frame(h_scores)
colnames(h_scores) <- c('h1', 'h2', 'h3')
h_scores <- cbind(wing_data_wide[,c(1:7)], h_scores)

# write.csv(h_scores, file =paste(".",outdir_tab,"h_scores.csv", sep="/"),
# row.names = FALSE, quote = TRUE)


# Section 4. Eigentensor Analysis --------------------------
# Tensor Analysis for REML estimates ONLY -
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
  REML_emat <- matrix(0,n,n)
  lowerTriangle(REML_emat) <- 1/sqrt(2) * REML_S_eigvec[(n+1):neigten,i]
  REML_emat <- REML_emat + t(REML_emat)
  diag(REML_emat) <- REML_S_eigvec[1:n,i]
  REML_S_eigTenmat[,,i] <- REML_emat
}

REML_S_eigten<-data.frame(NULL)
for (i in 1:neigten){REML_S_eigten<-rbind(REML_S_eigten,REML_S_eigTenmat[,,i])}

REML_S_eigtenvecs<-matrix(nrow=n*neigten, ncol=n)
REML_S_eigtenvals<-matrix(nrow=n*neigten, ncol=1)
for (i in 1:neigten){ #Using Emma's Method
  REML_S_eigtenvecs[((i-1)*n+1):(i*n),]=t(eigen(REML_S_eigten[((i-1)*n+1):(i*n),])$vectors) # eigenvectors in rows!!!!!!
  REML_S_eigtenvals[((i-1)*n+1):(i*n),1]=eigen(REML_S_eigten[((i-1)*n+1):(i*n),])$values
}


# Let's look at the Coordinates of the M matrix 
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

# Section 5. Eigentensor REML-MVN sampling -------------------
# Create Dave's Fig. 4 
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


MVN_S_val <- matrix(0, MVNsample, neigten)
colnames(MVN_S_val) <- paste("E", 1:neigten, sep="")
for (i in 1:MVNsample){
  for(j in 1:neigten){
    MVN_S_val[i,j] <- t(REML_S_vec[,j]) %*% MVN_S[,,i] %*% REML_S_vec[,j]
  }
}


# Get the null distribution for the eigentensor

n_cores <- detectCores()-1
cluster <- makeCluster(n_cores)
doParallel::registerDoParallel(cluster)
acomb <- function(...) abind(..., along=3)


null_eigTen_array <- foreach(i = 1:nullnumber, 
                        .packages = c('gdata'), 
                        .combine='acomb', .multicombine=TRUE) %dopar% {
                          # create null S matrix
                          temp_null_S<-array(NA, dim=c(neigten,neigten))
                          
                          # call in null estimates                   
                          temp_null_list <- asplit(null_M_array[,,i,], 3)
                          temp_null_list <- Filter(function(a) any(!is.na(a)), temp_null_list)
                          null_Mvarmat <- t(sapply(temp_null_list, diag))
                          null_Mcovmat <- t(sapply(temp_null_list, lowerTriangle))
                          
                          # fill out the S matrix
                          temp_null_S[1:n, 1:n]<- cov(null_Mvarmat,null_Mvarmat) #upper left quarter of S
                          temp_null_S[(n+1):neigten, (n+1):neigten]<-2*cov(null_Mcovmat,null_Mcovmat)#lower right quarter of S
                          temp_null_S[1:n, (n+1):neigten]<-sqrt(2)*cov(null_Mvarmat, null_Mcovmat)
                          temp_null_S[(n+1):neigten, 1:n]<-sqrt(2)*cov(null_Mcovmat, null_Mvarmat)
                          temp_null_S
                        }

parallel::stopCluster(cluster)

dim(null_eigTen_array)

# Get the average null S matrix and calcualte the eigenvectors
avg_null_S <- apply(null_eigTen_array, 1:2, mean)
avg_null_S_vecs <- eigen(avg_null_S)$vectors

null_S_eigvals <- laply(asplit(null_eigTen_array, 3), 
                        function(mat) diag(t(avg_null_S_vecs) %*% mat %*% avg_null_S_vecs))
null_S_CI <- t(apply(null_S_eigvals, 2, function(x) rangeFunc90(x)[2:3]))

null_S_CI <- cbind(as.integer(row.names(null_S_CI)),
                   eigen(avg_null_S)$values, null_S_CI)
colnames(null_S_CI)<- c("Eigten_Num","Var", "lowCI", "upCI")
null_S_CI <- as.data.table(null_S_CI)
null_S_CI$Eigten = paste0("E",null_S_CI$Eigten_Num)
null_S_CI$model <- "Null"



# Plot the 90% CIs for eigentensors and the NULL
#Add REML-MVN CI
EigtenCIs<-as.data.frame(NULL)
for (i in 1:neigten){
  EigtenCIs<-rbind(EigtenCIs, c(i, rangeFunc90(MVN_S_val[,i])))
}
colnames(EigtenCIs)<-c("Eigten_Num", "n", "lowCI", "upCI")
EigtenCIs %<>% mutate(Eigten=paste0("E", EigtenCIs$Eigten_Num))

#add reml estimates
REML_S_eigval %>% as.data.frame %>% 
  magrittr::set_colnames("REML_Var") %>%  
  mutate(Eigten_Num=row.names(.)) %>% 
  mutate(Eigten_Num=as.numeric(Eigten_Num)) %>%
  left_join(EigtenCIs, by=c("Eigten_Num"))->plot_EigtenCIs

# Add null estimates
plot_EigtenCIs <- as.data.table(plot_EigtenCIs)
plot_EigtenCIs$model <- "Obs"
colnames(plot_EigtenCIs)[1] <- "Var"

plot_EigtenCIs <- rbind(plot_EigtenCIs[,c("model","Eigten", "Eigten_Num", "Var",  "lowCI","upCI")],
                        null_S_CI[,c("model","Eigten", "Eigten_Num", "Var",  "lowCI","upCI")])
plot_EigtenCIs[, x_axis:=fcase(model == "Null", Eigten_Num + 0.15,
                               model == "Obs", Eigten_Num - 0.15)]

# Plot the variances explained by the Eigentensors with the REML-MVN CIs
# cairo_ps(paste(".",outdir_fig,"Eigentenor_Var.eps", sep="/"),
#            family = "Times", pointsize=14, width =6, height = 4.35)

# png(paste(".",outdir_fig,"Eigentenor_Var.png", sep="/"),  res=300,
#     bg = "white", pointsize=12, width = 1340.5, height = 1126.5)

par(mfrow = c(1,1))
Eigtenlim = 11
par(mar = c(2.5,4.5,0.5,0.5))
plot(NULL, type = "n", xlab = "", yaxs = "i", ylab = "Among-Line Variance", 
     main = "",  xaxt = "n", yaxt = "n",
     xlim = c(0.5,Eigtenlim+0.5), ylim = c(0,0.06), bty = "L")
epsilon=0.1
for (mod in c("Obs", "Null")) {
  plot_EigtenCIs[model== mod& Eigten_Num <= Eigtenlim, 
                 segments(x_axis, lowCI, x_axis, upCI, 
                          lty = ifelse(mod == "Obs",1,2)) ]
  plot_EigtenCIs[model== mod & Eigten_Num <= Eigtenlim, 
                 segments(x_axis-epsilon, lowCI, x_axis+epsilon, lowCI)]
  plot_EigtenCIs[model== mod & Eigten_Num <= Eigtenlim, 
                 segments(x_axis-epsilon, upCI, x_axis+epsilon, upCI) ]
  plot_EigtenCIs[model== mod & Eigten_Num <= Eigtenlim,
                 points(x_axis, Var, pch =  ifelse(mod == "Obs",16,21))]
}
axis(2, at=seq(0,0.06, 0.01), labels=sprintf("%.2f", seq(0,0.06, 0.01)), las=2)
axis(1, at=c(1:Eigtenlim), 
     labels=as.expression(lapply(1:Eigtenlim, function(i)bquote(bold("E")[.(i)]))))
legend("topright", c("Observed", "Randomised"), bty="n", 
       pch=c(16,21), lty=c(1,2))

# dev.off()



# Get Major axes of Eigentensors
Project_major_EigtenVecs <- foreach(evec = c(1,6,7,13,18),.combine=rbind) %do% {
  e_vec <- REML_S_eigtenvecs[evec,]
  enum = ifelse(evec<7, (10+evec), ifelse(evec>10, evec+18, 21))
  # Get CIs for M
  e_vec_thrmat <- t(apply(apply(AsycovM_array, c(3,4), 
                                function (x) t(e_vec) %*% x %*% e_vec),2,function (t) rangeFunc90(t)))
  e_vec_thrmat <- cbind(rep(enum,12), apply(M_array, 3, 
                                            function (x) t(e_vec) %*% x %*% e_vec), e_vec_thrmat)
  
}; rm(e_vec_thrmat)
colnames(Project_major_EigtenVecs) <- c("vec", "var", "n", "loCI", "upCI")
matInfo <- rownames(Project_major_EigtenVecs)
Project_major_EigtenVecs <- as.data.frame(Project_major_EigtenVecs)
setDT(Project_major_EigtenVecs)
Project_major_EigtenVecs$matInfo <- matInfo
Project_major_EigtenVecs[, c("Treat","Gen") := data.table(str_split_fixed(matInfo,"_",2))]
Project_major_EigtenVecs[, x:=fcase(Treat=="1", as.numeric(Gen)-0.15,
                                    Treat=="3", as.numeric(Gen)+0.15,default=NA)]

# cairo_ps(paste(".",outdir_fig,"Eigenten_eigvec_proj.eps", sep="/"),
#            family = "Times", pointsize=11, width =6, height = 4.35)

png(paste(".",outdir_fig,"Eigenten_eigvec_proj.png", sep="/"), res=300,
    bg = "white", pointsize=14, width = 2044, height = 1441.533)

par(mfrow = c(2,3), mar=c(3,3,2.5,1), oma=c(2,2,0,0))
for (evec in c(11,16,21,31,36)) {
  Project_major_EigtenVecs[vec ==  evec,
                           plot(x,var, xlim=c(0.5, 6.5), xlab="",las=2,
                                ylim = c(min(loCI), max(upCI)),xaxt="n",
                                main=bquote(italic("e")[.(evec)]),
                                ylab="")]
  abline(h=0, lty=2)
  axis(side=1, at=1:6, labels = 1:6)
  Project_major_EigtenVecs[vec == evec,segments(x,loCI,x, upCI)]
  Project_major_EigtenVecs[vec == evec,segments(x-0.075,loCI,x+0.075, loCI)]
  Project_major_EigtenVecs[vec == evec,segments(x-0.075,upCI,x+0.075, upCI)]
  Project_major_EigtenVecs[vec == evec,
                           points(x,var, pch=21, bg=ifelse(Treat=="1","black", "white"))]
  
}
mtext("Among-line varaince \u00B1 90% CI", side=2, line= 0.75, outer = T, cex=0.85)
mtext("Generations", side=1, at = 3.5, line= 3.25, cex=0.85)
plot.new()
legend("center", pch =c(16, 21), col = c(1,1),
       bg =c(1,0),
       legend= c("Small", "Large"), bty="n")
dev.off()


#Table for Supplementary
as.data.frame(apply(matrix(REML_S_eigtenvals[1:(11*6),], ncol=6,  byrow =TRUE), 
      2, function(x) sprintf("%1.3f", x)))->REML_S_eigtenvals_Tab
REML_S_eigtenvals_Tab<-cbind(c(1:11),REML_S_eigtenvals_Tab)
colnames(REML_S_eigtenvals_Tab)<-c("Ek", paste0("e",1:6))
# write.table(REML_S_eigtenvals_Tab, file ="C:/Users/carac/Dropbox/Analysis/Multivariate_MR/Tables/eigentensor_eigenvalues.txt",
#           row.names = FALSE, quote = TRUE, sep="\t")


# r Vector correlation
# vectors flipped so major loading within eigenvector is positive
mjrEvecs <- t(REML_S_eigtenvecs[c(1,6,7,13,18),]) %*% diag(c(-1,-1,1,-1,1)) 

# Dotproduct table
DotProd<-rbind(apply(mjrEvecs, 2, function (x) t(H_vecs[,1])%*% x),
               apply(mjrEvecs, 2, function (x) t(H_vecs[,2])%*% x),
               apply(mjrEvecs, 2, function (x) t(H_vecs[,3])%*% x))
apply(DotProd,2, function (x) round(x,digits=2)) %>% View()

# apply(DotProd,2, function (x) rad2deg(acos(x)))


# Section 6. Coordinates analyses ######
# get confidence intervals on the VL Coordinates using the 10,00 REML-MVN estimates of M
# set up Coords array
VLcoordsMVN <- array(NA,dim=c(1,1,p,neigten,MVNsample))

for (i in 1:MVNsample) {
  for (j in 1:neigten) {
    for (k in 1:p){
      frobenius.prod(REML_S_eigTenmat[,,j], AsycovM_array[,,i,k])->VLcoordsMVN[,,k,j,i]
    }
  }
}

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

setDT(CoordVL_CIs_fin)

# 
# png(paste(".",outdir_fig,"VL_Coords.png", sep="/"), res=300,
#     bg = "white", pointsize = 14, width = 797.25, height = 2162.3)

par(mfrow = c(3,1), mar = c(3,3,2,1), oma = c(2,2,0,0))
for (i in c("C_E1", "C_E2", "C_E3")) {
  C_num = gsub("C_E", "", unique(CoordVL_CIs_fin[cord_Eig == i]$cord_Eig))
  # We flipped E2
  if (i != "C_E2"){
    Cord_dat <- CoordVL_CIs_fin[cord_Eig == i]
    Cord_dat[, plot(xnum_coords, Coord, las=2, xaxt="n",
                    main = bquote(bold(E)[{.(C_num)}]),pch=16,
                    xlab="", xlim=c(0.5,6.5), 
                    ylim = c(floor((min(lowCI)*5))/5, ceiling((max(upCI)*5))/5),
    )]
  } else {
    Cord_dat <- CoordVL_CIs_fin[cord_Eig == i, .(Coord,lowCI, upCI)]
    Cord_dat <- cbind(CoordVL_CIs_fin[cord_Eig == i,.(pop, xnum_coords)],(Cord_dat * -1))
    Cord_dat[, plot(xnum_coords, Coord , las=2, xaxt="n",pch=16,
                        main = bquote(bold(E)[{.(C_num)}]),
                        xlab="", xlim=c(0.5,6.5), 
                        ylim = c(-0.2, 0.6),
    )]
  }
  abline(h=0, lty=2)
  axis(side=1, at = 1:6, labels = 1:6)
  Cord_dat[, segments(xnum_coords-0.1, lowCI, xnum_coords+0.1, lowCI)]
  Cord_dat[, segments(xnum_coords-0.1, upCI, xnum_coords+0.1, upCI)]
  Cord_dat[, segments(xnum_coords, lowCI, xnum_coords, upCI)]
  Cord_dat[!pop %in% c(1:6), points(xnum_coords, Coord, pch=21, bg="white")]
  if (i == "C_E1"){
    legend("topright", pch =c(16, 21), bg="white", legend = c("S", "B"), bty = "n")
  }
}
mtext("Coordinates \u00B1 90% CI", side = 2, cex = 0.85, line = 0.25, outer = T )
mtext("Generations", side = 1, at=3.5, cex = 0.85, line = 3.25 )

# dev.off()

# Section 7. Variability as CVs  ------------------
# Get first three eigenvalues of M
trait_CVs <- unlist(lapply(M_list, function(x) eigen(x)$values[1:3]))
trait_CVs <- as.data.frame(trait_CVs)
trait_CVs$pop <- row.names(trait_CVs)
trait_CVs$trait <- paste0("m", substr(trait_CVs$pop, 4, 4))
trait_CVs$pop <- substr(trait_CVs$pop, 1, 3)
names(trait_CVs)[1] <- "vl"

# Get M traces
trace_CV <- as.data.frame(unlist(lapply(M_list, tr)))
trace_CV <- as.data.frame(trace_CV)
trace_CV$pop <- row.names(trace_CV)
trace_CV$trait <- "trace"
names(trace_CV)[1] <- "vl"

trait_CVs <- rbind(trait_CVs, trace_CV); rm(trace_CV)

setDT(trait_CVs)


# Add in h vl
h_vl_temp <- data.frame("vl"= H_thru_M$var, 
                        "pop" = H_thru_M$pop, 
                        "trait" = paste0("h", H_thru_M$Hvec)) 

trait_CVs <- rbind(trait_CVs,
                   h_vl_temp); rm(h_vl_temp)

# add in e vl
e_vl_temp <- data.frame("vl"= Project_major_EigtenVecs$var, 
                        "pop" = Project_major_EigtenVecs$matInfo, 
                        "trait" = paste0("e", Project_major_EigtenVecs$vec)) 

trait_CVs <- rbind(trait_CVs,
                   e_vl_temp); rm(e_vl_temp)


# calculate CVs
trait_CVs <- trait_CVs[, (sd(vl)/mean(vl))*100, by="trait"]
names(trait_CVs)[2] <- "CV"

# Add the trait variances, calculating CV per trait
tr_var_temp <- M_tab_lng[trait2 == trait, .(vl, pop, trait)]

trait_CVs <- rbind(trait_CVs,
                   data.frame("trait" = "var", 
                              "CV" = tr_var_temp[, (sd(vl)/mean(vl))*100, by="trait"]$V1))
# Tidy up names
trait_CVs$vec <- gsub('[0-9]+', '', trait_CVs$trait)
trait_CVs$num <- gsub(".*?([0-9]+).*", "\\1", trait_CVs$trait)
trait_CVs[trait %in% c("trace", "var"), c("vec", "num"):= ""]
trait_CVs[, x:= .I]
trait_CVs[trait == "var", x:= 0]

dev.off()

# png(paste(".",outdir_fig,"trait_est_cv.png", sep="/"),  res=300,
#     bg = "white", pointsize=12, width = 2244, height = 987)

# plot
par(mfrow = c(1,1), mar=c(4,4,0,0.5), oma= c(0,0,0,0))
plot(NULL, xlim=c(-0.5,12.5), xlab="Trait relationship", bty="L",
     ylab="Coefficient of variance",xaxt="n",yaxs="i", ylim=c(10,90), las = 2)
trait_CVs[trait == "var",
          boxplot(CV~ trait, col = "white",  
                  yaxt = "n", xaxt = "n", at = 0,
                  frame = F,boxwex = 1, add = TRUE)]
for (i in 1:12) {
  # trait_CVs[x==i, boxplot(CV ~ trait, col = "white",  
  #                   yaxt = "n", at = i,  xaxt = "n",
  #                   frame = F,boxwex = 1, add = TRUE)]
  trait_CVs[x == i,points(x, CV, pch=16)]
  trait_CVs[x == i, axis(side = 1, at = i, 
                       labels = bquote(italic(.(vec))[.(num)]))]
}
# points(0, mean(trait_CVs[trait=="var"]$CV), pch=16)
axis(side=1, at = c(0,4), labels = c("variance", "trace"))
# dev.off()
















# 
# # Section 8. Mitteroecker Distance ######
# # pr.coord (principle co-ordinates analysis)
# # doesn't work for R 
# # doesn't work for Eigvectors as one M is not positive definite
# library(vcvComp)
# citation("vcvComp")
# 
# Hspace_arrayM<-array(NA, dim=c(3,3,12))
# for(i in 1:12){
#   Hspace_arrayM[,,i]<-t(H_vecs) %*% M_array[,,i] %*% H_vecs
# }
# 
# # check matrices are pos.def
# apply(Hspace_arrayM[,,],3,is.positive.definite)
# 
# 
# # Now, to calculate the D (distance matrix) and get the principle 
# disMAtH_M<-matrix(NA, nrow = 12, ncol = 12)
# for (i in 1:12){
#   k=i
#   while (k<=12) {
#     Hspace_arrayM[,,i] %*% solve(Hspace_arrayM[,,k])->tmp
#     sqrt(sum(log(eigen(tmp)$values)^2))->disMAtH_M[i,k]
#     k=k+1
#   }
# }
# lowerTriangle(disMAtH_M) = upperTriangle(disMAtH_M, byrow=TRUE)
# 
# # pr.coord (principle co-ordinates analysis)
# vcvComp::pr.coord(disMAtH_M)->prcoordSupH_M
# 
# 
# # Scree plot of principle coordinates
# plot_prcoordSupH<-as.data.frame(cbind(1:11,prcoordSupH_M$values[1:11]))
# colnames(plot_prcoordSupH)<-c("Cord","Eigval")
# setDT(plot_prcoordSupH)
# plot_prcoordSupH[,Eigval:=as.numeric(Eigval)]
# plot_prcoordSupH[,Prop:=round(Eigval/sum(Eigval)*100, 1)]
# 
# as.data.frame(cbind(1:12, prcoordSupH_M$PCoords[,1:3]))->PlotPcords
# colnames(PlotPcords)[1]<-c("p")
# PlotPcords$pop=unlist(names(M_list))
# setDT(PlotPcords)
# PlotPcords[, c("Treat","Gen") := data.table(str_split_fixed(pop,"_",2))]
# PlotPcords[, c("Treat","Gen") := PlotPcords[, lapply(.SD, as.numeric), .SDcols =c("Treat","Gen")]]
# 
# #plot
# # svglite(filename = "C:/Users/carac/Dropbox/Analysis/Multivariate_MR/Multi_Latex_Documents/PNGs and EPS/M_in_H_MittDist.svg",
# #     width = 7.5, height = 5.4, pointsize = 11,
# #     bg = "white", system_fonts = "Times New Roman")
# # grDevices::cairo_pdf(filename =
# #                       "C:/Users/carac/Dropbox/Analysis/Multivariate_MR/Multi_Latex_Documents/PNGs and EPS/M_in_H_MittDist.pdf",fallback_resolution = 2400,
# #                     family = "Times", bg="transparent",pointsize=11, width = 7.5, height = 5.4)
# 
# #including a 3d plot
# # grDevices::cairo_pdf(filename =
# #                       "C:/Users/carac/Dropbox/Analysis/Multivariate_MR/Multi_Latex_Documents/PNGs and EPS/M_in_H_MittDist_3D.pdf",fallback_resolution = 2400,
# #                     family = "Times", bg="transparent",pointsize=11, width = 7.25, height = 8.4)
# # svglite(filename = "C:/Users/carac/Dropbox/Analysis/Multivariate_MR/Multi_Latex_Documents/PNGs and EPS/M_in_H_MittDist_3D.svg",
# #     width = 7.25, height = 8.4, pointsize = 11,
# #     bg = "white", system_fonts = "Times New Roman")
# m<-matrix(c(rep(1:2, each = 3),
#   rep(1:2, each = 3),
#   rep(3:4, each = 3),
#   rep(3:4, each = 3),
#   rep(5, 18)), 7,6, byrow = TRUE)
# layout(m)
# par(mar = c(5,5,2,1))# mfrow=c(2,2),
# plot_prcoordSupH[,plot(Cord, Prop, bty="L",ylim=c(0,40),
#                        ylab="", xlab="Principal coordinate",
#                        xlim=c(0,12),xaxs="i",family = "Times New Roman",
#                        yaxs="i",pch=16, xaxt="n", yaxt="n")]
# axis(1, at=c(seq(1,11,1)), labels=seq(1,11,1),family = "Times New Roman", cex.axis=1)
# axis(2, at=c(seq(0,40, 5)), las=2,cex=0.6, family = "Times New Roman",
#      labels=paste0(format(round(seq(0,40, 5),1), nsmall = 0),"%"))
# axis(2, at=20, labels = "Proportion of varaince",family = "Times New Roman", line=2.1, tick=F)
# # plot_prcoordSupH[,text(Cord, Prop, Prop, pos=3, cex=0.8)]
# mtext("A", 3,family = "Times New Roman", outer=FALSE, cex=1.35,adj=-0.22, line=.65)
# # 47% of the information captured in the first 2 principle coordinates, and 58% in 3
# #  2d with lines
# par(xpd=TRUE, family = "Times New Roman")
# for (i in 1:3) {
#   if (i==1){
#     Cols<-c("PCo1","PCo2")
#     xlimy=c(-1,2);ylimy=c(-1,1)
#     lab<-c("First principal coordinate", "Second principal coordinate")
#     labpos_1=c(2,2,1,1,2,4)
#     labpos_3=c(4,2,1,1,4,2)
#   }else if (i==2){
#     Cols<-c("PCo1","PCo3")
#     xlimy=c(-1,2);ylimy=c(-1.5,1)
#     lab<-c("First principal coordinate", "Third principal coordinate")
#     labpos_1=c(4,2,4,4,3,4)
#     labpos_3=c(4,3,1,3,1,1)
#   }else {
#     Cols<-c("PCo2","PCo3")
#     xlimy=c(-1,1.5);ylimy=c(-1.5,1)
#     lab<-c("Second principal coordinate", "Third principal coordinate")
#     labpos_1=c(1,2,1,2,3,2)
#     labpos_3=c(4,2,2,1,3,4)
#   }
#   plot(NULL, ylim=ylimy, yaxt="n", ylab=lab[2],family = "Times New Roman",
#        xlab=lab[1], bty="L", xlim=xlimy)
#   axis(side=2, at=seq(ylimy[1],ylimy[2],0.5), las=2, family = "Times New Roman",
#        labels=sprintf("%1.1f",seq(ylimy[1],ylimy[2],0.5)))
#   abline(h=0, lty=2, xpd=F);abline(v=0, lty=2, xpd=F)
#   PlotPcords[Treat==1, lines(.SD, lwd=2, 
#                              col=rgb(0.2,0.2,0.2,alpha=0.75)),.SDcols = Cols]
#   PlotPcords[Treat==1, points(.SD, pch=16, cex=1.15),.SDcols = Cols]
#   PlotPcords[Treat==3, lines(.SD, lwd=2, col=rgb(0.1,0.2,1,alpha=0.75)),.SDcols = Cols]
#   PlotPcords[Treat==3, points(.SD, pch=16, cex=1.15, col="blue"),.SDcols = Cols]
#   PlotPcords[Treat==1,text(.SD, labels=Gen, family = "Times New Roman",pos=labpos_1),.SDcols = Cols]
#   PlotPcords[Treat==3,text(.SD, labels=Gen, family = "Times New Roman",pos=labpos_3, col="blue"),.SDcols = Cols]
#   mtext(LETTERS[i+1], 3,family = "Times New Roman", outer=FALSE, cex=1.35,adj=-0.22, line=.65)
# }
# legend("bottomright",pch=16, col=c(1,"blue"),lty=1,
#        c("Small population", "Large Population"),bty="n")
# # dev.off()
# 
# #Offset labels
# par(mar = c(2,5,2,1))
# PC_off<-PlotPcords[,.(p, pop, Treat, Gen)]
# PC_off$x=c(0,-0.15,0,-0.15,0,0.06, 0.08,-0.12,0,0,0.07,0)
# PC_off$z=c(-0.09,0,-0.11,0,0.09,0, 0,0,-0.09,-0.07,0,0.08)
# 
# library("scatterplot3d")
# PlotPcords_3D<-merge(PlotPcords, PC_off)
# PlotPcords_3D %>% as.data.frame %>% 
#   with(scatterplot3d(x=PCo1,y=PCo2,z=PCo3,
#                             type = "h",pch=16,
#                             xlab="First principal coordinate", 
#                             ylab="", color=ifelse(Treat==1, "black", "blue"),
#                             zlab="Third principal coordinate"))
# text(6.5, -2, "Second principal coordinate", srt=50)
# 
# library("plot3D")
# PlotPcords_3D<-merge(PlotPcords, PC_off)
# PlotPcords_3D[Treat==1,scatter3D(x=PCo1,y=PCo2,z=PCo3,pch=16,nticks =6, label=TRUE,
#                      col="black", type="h",colkey = FALSE, ticktype="detailed", 
#                      axes=TRUE,theta=30, phi=25,bty = "b2",
#                      xlim=c(-1, 1.8), ylim=c(-1, 1),zlim=c(-1.1, 0.6),
#                      xlab="First principal coordinate", 
#                      ylab="Second principal coordinate", 
#                      zlab="Third principal coordinate") ]
# PlotPcords_3D[Treat==1,scatter3D(x=PCo1,y=PCo2,z=PCo3,pch=16, label=TRUE,
#                      col="black",add = TRUE, type="b", lwd=2)]
# PlotPcords_3D[Treat==3,scatter3D(x=PCo1,y=PCo2,z=PCo3,pch=16, label=TRUE,
#                      col="blue", add = TRUE, type="h")]
# PlotPcords_3D[Treat==3,scatter3D(x=PCo1,y=PCo2,z=PCo3,pch=16, label=TRUE,
#                      col="blue",add = TRUE, type="b", lwd=2)]
# PlotPcords_3D[,text3D(x=PCo1+x,y=PCo2,z=PCo3+z, colkey = FALSE, add = TRUE, 
#         labels = Gen, colvar=Treat, col = c("black", "blue"))]
# legend(0.5,0,pch=16, col=c(1,"blue"),lty=1,
#        c("Small population", "Large Population"),bty="n")
# mtext(LETTERS[5], side=3, cex=1.35, adj=0.25,family = "Times New Roman",line=-0.35)
# # dev.off()

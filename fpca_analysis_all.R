####################################################################################
### Functional principle component analysis (fPCA) on Actigraphy accelerometer data
### R code prepared for manuscript: 
####################################################################################
#### R code accelerometer data from NHANES 2011-2014 and perform the following
#### 1. calculate the mean activity value at each 5-min epoch across days for each subject
#### 2. perform fPCA analysis
#### 3. fit functional object split by median of fPCA eigenvalues

### load libraries
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(dplyr); library(fda);library(tidyr);library(mgcv); library(purrr) 

### load script for different color scheme for making figures
source("plot_fd_new1.r")
source("plot_fpca_new1.r")
####################################################################################
#### 1. calculate the mean activity value at each 1-min epoch across days for each subject
####################################################################################
#### read in sample data
data0 <- readRDS("/Users/dsuolang/Desktop/Study1/Study1/nhanes_acc_clean2_10min.rds")
#data0$tenmin<-data0$fivemin
#data0$fivemin<-NULL
ID <- unique(data0$SEQN)
dat_wide <- data0[order(data0$SEQN), ]
dat_wide <- pivot_wider(data0, names_from = SEQN, values_from = activity_avg)
dat_wide$day <- NULL
new_colnames <- c("tenmin", paste0("ID", names(dat_wide)[-1])) #change
colnames(dat_wide) <- new_colnames
# Split the data into subsets of size 1440 and add suffix to column names
split_data <- split(dat_wide, rep(1:7, each = 144))
combined_data <- do.call(cbind, split_data)
combined_data <- combined_data[, -grep("tenmin", colnames(combined_data))]
combined_data$tenmin<-1:144

dat_wide<-NULL
data0<-NULL
####################################################################################
#### 2. perform fPCA analysis
####################################################################################
dat <- combined_data
ID <- substr(colnames(dat)[-1],1,11)
length(ID)
### setup fitting functional object
n_basis <- 9 # number of basis functions, can only be odd numbers
### timepoints for x axis
timepoints <- 1:144 #change
## create fourier basis functions
basis <- create.fourier.basis(c(1,length(timepoints)), nbasis=n_basis)
## fit functional object
fd <- smooth.basis(timepoints, data.matrix(dat[,-1]), basis, method = "qr")$fd
L <- dim(fd$y)[1] # get the number of data samples.
fpar <- fdPar(basis) # functional parameter object
### first fit fPCA -- preliminary, to decide how many components to use
### here to we extract the first 9 PCA and then examine the scree plot
fpca0 <- pca.fd(fd, nharm=9, fpar)
# Calculate margins based on the number of digits, set margins
sig_digits <- max(floor(log10(fpca0$values)), na.rm = TRUE)
left_margin <- 1 + (sig_digits + 1) * 1 # Adjust these values as needed
bottom_margin <- 3
top_margin <- 3
right_margin <- 1
par(mar = c(bottom_margin, left_margin, top_margin, right_margin))

plot(fpca0$values, type = "b", ylab = "eigenvalue", xlab = "number of components")
## choose an appropriate number of components to use in the final analysis
### here we choose 4 PCA based on the scree plot
fpca <- pca.fd(fd,nharm=4, fpar)
### plot fPCA components
## mean activity is same in all figures
## "+" shows add the component to the mean
## "-" shows subtract the component to the mean
par(mfrow=c(2,2))
plot(fpca,cex.main=.8)
flip <- c(1,2,4) # identify which components to flip the sign
plot.fpca.new(fpca,cex.main=0.8, flip = flip) # flipped fPCA
### extract corresponding eigenvalues
score0 <- fpca$scores
# flip the sign for fPCA1,2,4
score0[,flip] <- -score0[,flip]
fpca$scores <- score0
# DO not standardize with mean and sd
#score <- scale(fpca$scores)
score1 <- data.frame(ID = ID, score = fpca$scores)
colnames(score1)[2:5] <- paste0("fPCA", 1:4)
# Save the results for the current day
#day_result_filename <- paste0("score_df_all_10min_2024.rds") #change
#saveRDS(score1, file = day_result_filename)

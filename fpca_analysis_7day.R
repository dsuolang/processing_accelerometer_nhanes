####################################################################################
### Functional principle component analysis (fPCA) on Actigraphy accelerometer data
### R code prepared for manuscript: 
####################################################################################
#### R code used 50 randomly selected subjects with accelerometer data from NHANES 2011-2014 and perform the following
#### 1. calculate the mean activity value at each 5-min epoch across days for each subject
#### 2. perform fPCA analysis
#### 3. fit functional object split by median of fPCA eigenvalues

### load libraries
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(dplyr); library(fda);library(tidyr);library(mgcv) 

### load script for different color scheme for making figures
source("plot_fd_new1.r")
source("plot_fpca_new1.r")

####################################################################################
#### 1. calculate the mean activity value at each 5-min epoch across days for each subject
####################################################################################
#### read in sample data
data0 <- readRDS("/Users/dsuolang/Desktop/Study1/Study1/nhanes_acc_clean2_10min.rds")# actual full data
#data0$tenmin<-data0$fivemin
#data0$fivemin<-NULL
ID <- unique(data0$SEQN)
# Specify the number of days you want to process
num_days <- 7

# Loop over the days
for (day in 1:num_days) {
  # Select rows where 'day' is equal to the current day
  day_data <- data0[data0$day == day, ]
  #reshape
  dat_wide <- pivot_wider(day_data, names_from = SEQN, values_from = activity_avg)
  dat_wide$day <- NULL
  new_colnames <- c("tenmin", paste0("ID", names(dat_wide)[-1]))
  colnames(dat_wide) <- new_colnames
  dat_wide <- head(dat_wide, 144)
  
  #dat_wide<-NULL
  #data0<-NULL
  ####################################################################################
  #### 2. perform fPCA analysis
  ####################################################################################
  ID <- substr(colnames(dat_wide)[-1], 3, 8)
  ### setup fitting functional object
  n_basis <- 9 # number of basis functions, can only be odd numbers
  ### timepoints for x axis
  timepoints <- 1:144 #change
  ## create fourier basis functions
  basis <- create.fourier.basis(c(1,length(timepoints)), nbasis=n_basis)
  ## fit functional object
  fd <- smooth.basis(timepoints, data.matrix(dat_wide[,-1]), basis, method = "qr")$fd
  L <- dim(fd$y)[1] # get the number of data samples.
  fpar <- fdPar(basis) # functional parameter object
  ### first fit fPCA -- preliminary, to decide how many components to use
  ### here to we extract the first 9 PCA and then examine the scree plot
  fpca0 <- pca.fd(fd, nharm=9, fpar)
  # Calculate margins based on the number of digits, set margins
  sig_digits <- max(floor(log10(fpca0$values)), na.rm = TRUE)
  left_margin <- 1 + (sig_digits + 1) * 1 # Adjust these values as needed
  bottom_margin <- 2
  top_margin <- 2
  right_margin <- 1
  par(mar = c(bottom_margin, left_margin, top_margin, right_margin))
  #plot(fpca0$values, type = "b", ylab = "eigenvalue", xlab = "number of components")
  ## choose an appropriate number of components to use in the final analysis
  ### here we choose 4 PCA based on the scree plot
  fpca <- pca.fd(fd,nharm=4, fpar)
  ### plot fPCA components
  ## mean activity is same in all figures
  ## "+" shows add the component to the mean
  ## "-" shows subtract the component to the mean
  par(mfrow = c(1, 4), pin = c(1.15, 1.15))
  #plot(fpca,cex.main=.8)
  flip <- c(1,2,4) # identify which components to flip the sign
  plot.fpca.new(fpca,cex.main=.8, flip = flip) # flipped fPCA
  ### extract corresponding eigenvalues
  score0 <- fpca$scores
  # flip the sign for fPCA1,2,4
  score0[,flip] <- -score0[,flip]
  fpca$scores <- score0
  # Standardize with mean and sd
  score <- scale(fpca$scores)
  score1 <- data.frame(ID = ID, score = score)
  colnames(score1)[2:5] <- paste0("fPCA", 1:4)
  # Save the results for the current day
  day_result_filename <- paste0("score1_df_day", day, ".rds")
  saveRDS(score1, file = day_result_filename)
}

# DO NOT RUN BELOW


# Define the directory where your RDS files are located and where you want to save the merged data
# data_dir <- "home/dsuolang/Study1"  # Replace with your actual directory path
# # Specify the number of days you want to process
# num_days <- 7
# # Create an empty list to store the data frames
# result_list <- list()
# # Load and merge data from individual files
# for (day in 1:num_days) {
#   day_result_filename <- file.path(data_dir, paste0("score1_df_day", day, ".rds"))
#   day_data <- readRDS(day_result_filename)
#   
#   result_list[[day]] <- day_data
# }
# # Combine all data frames into one based on the 'ID' column
# merged_fpca <- do.call(cbind, result_list)
# merged_fpca<- merged_fpca[ , !duplicated(names(merged_fpca))]
# # Assuming 'ID' is a variable in the 'merged_fpca' data frame
# merged_fpca <- merged_fpca %>%
#   rename(seqn = ID)
# saveRDS(merged_fpca,"merged_fpca.rds")

### set directory
# setwd("your directory")
setwd("D:/log/대학원/연구/졸업논문/Data")

### load library
# cox regression
library(fda)
library(survival)
library(grpreg)

# data preprocessing
library(dplyr)

# brain Visualization
library(ggseg)
library(patchwork)
library(ggseg3d)
library(ggplot2)
library(rgl)
library(misc3d)
library(neurobase)
if (!requireNamespace("aal")) {
  devtools::install_github("muschellij2/aal")
  library(aal)
} else {
  library(aal)
}
if (!requireNamespace("MNITemplate")) {
  devtools::install_github("jfortin1/MNITemplate")
  library(MNITemplate)
} else {
  library(MNITemplate)
}

### set IBS kernel 
# X is matrix by Gene x Brain data 
# genenum is # of Gene data in X
# brainnum is # of Brain data in X

ibs.kernel.intigrate.nobrain <- function(X){
  tmparray <- matrix(0, nrow = nrow(X), ncol = nrow(X))
  for(j in 1:ncol(X)){
    N = nrow(X)
    colm <- matrix(rep(X[,j],N),nrow = N)
    rowm <- matrix(rep(X[,j],N),nrow = N,byrow = T)
    
    tmpmat <- ifelse(colm == rowm, 2, ifelse(abs(colm-rowm) == 1, 1, 0))/2
    tmparray <- tmparray + tmpmat
  }
  return(tmparray/j)
}

kernel.matrix.inti.nobrain <- function(X, generange, brainrange, method = "gaussian", rho = 1){
  ls = list()
  if(method == "gaussian"){
    for(j in generange){
      print(paste0("num ",j," in gene"))
      for(k in brainrange){
        print(paste0("num ",k," in roi"))
        assign(paste0("ker",j,k),gaussian.kernel.data(X, j, k, rho))
        ls <- mget(ls(pattern = "k"))
      }
    }
  }else if(method == "ibs"){
    genedata <- X[,generange]
    name <- names(X[,generange])
    gene <- unlist(strsplit(name,"_"))[seq(1,length(unlist(strsplit(name,"_"))),2)]
    for(i in unique(gene)){
      print(paste0(i," in gene"))
      onegenedata <- genedata[,which(gene==i)]
      assign(paste0("geker",i),ibs.kernel.intigrate.nobrain(onegenedata))
      ls <- mget(ls(pattern = "k"))
    }
  }else{
    cat("it has no method. You chagne another method. /n")
    return()
  }
  return(ls)
}

################################################################################################################
kernel.matrix.inti.withbrain <- function(X, generange, brainrange, method = "gaussian", rho = 1){
  ls = list()
  if(method == "gaussian"){
    for(j in generange){
      print(paste0("num ",j," in gene"))
      for(k in brainrange){
        print(paste0("num ",k," in roi"))
        assign(paste0("ker",j,k),gaussian.kernel.data(X, j, k, rho))
        ls <- mget(ls(pattern = "k"))
      }
    }
  }else if(method == "ibs"){
    genedata <- X[,generange]
    name <- names(X[,generange])
    gene <- unlist(strsplit(name,"_"))[seq(1,length(unlist(strsplit(name,"_"))),2)]
    for(i in unique(gene)){
      print(paste0(i," in gene"))
      onegenedata <- genedata[,which(gene==i)]
      for(j in brainrange){
        print(paste0("num ",j," in roi"))
        tptp = diag(X[,j]) %*% ibs.kernel.intigrate.nobrain(onegenedata) %*% diag(X[,j])
        assign(paste0("ker",i,j),tptp)
        ls <- mget(ls(pattern = "k"))
      }
    }
  }else{
    cat("it has no method. You chagne another method. /n")
    return()
  }
  return(ls)
}

### transfer list to dataframe
listfordataframe <- function(listk){
  kkk <- data.frame(0)
  for(i in 1:length(listk)){
    df <- as.data.frame(listk[[i]])
    kkk <- cbind(kkk, df)
  }
  kk <- kkk[,2:length(kkk)]
  names(kk) <- names(as.data.frame(listkernel))
  return(kk)
}

### transfer list to dataframe
listfordataframe.nobrain <- function(listk){
  kkk <- data.frame(0)
  for(i in 1:length(listk)){
    df <- as.data.frame(listk[[i]])
    kkk <- cbind(kkk, df)
  }
  kk <- kkk[,2:length(kkk)]
  names(kk) <- names(as.data.frame(test))
  return(kk)
}

#########################################################################################################################
### example my data
b <- read.csv("transdata.csv",header=T)

# checked this data column num in
## y 2, 3
## demographic 4:sex 5:age
## gene 6:25
## roi 26:59
listkernel <- kernel.matrix.inti.withbrain(b, 6:25, 26:59, method = "ibs") # 59
interaction.data <- listfordataframe(listkernel)

### fitting grpsurv
# test for gene data
test <- kernel.matrix.inti.nobrain(b, 6:25, method="ibs")
gene.kernel.data <- listfordataframe.nobrain(test)

# cbind of X
X <- b[,c(1:5,26:59)] %>%
  dplyr::select(-PATNO, -followupmon, -fog) %>%
  cbind(gene.kernel.data) %>% 
  cbind(interaction.data)
  
# make y
y <- b %>%
  dplyr::select(followupmon, fog)

# duplicate groups
group <- names(X)
group[1:36] <- "0" # no kernel

int.kernel.groups <- substr(names(interaction.data),1,9) # gene kernel
gene.kernel.groups <- substr(names(gene.kernel.data),1,9) # interaction kernel

group <- c(group[1:36], gene.kernel.groups, int.kernel.groups)

#################################################
### check this mulitplier
# SVD for Kernel
k <- c()
for(ii in 1:4){ # gene counts SNCA, GBA, LRRK2, COMT 4 
  # std for kernel matrix
  std_test <- .Call("standardize", as.matrix(test[[ii]]))
  
  # find # of eigen value
  k <- c(k, length(which(svd(std_test[[1]])$d > 1e-10)))
}

for(jj in 1:136){ # interaction counts 4*34
  # std for kernel matrix
  std_test <- .Call("standardize", as.matrix(listkernel[[jj]]))
  
  # find # of eigen value
  k <- c(k, length(which(svd(std_test[[1]])$d > 1e-10)))
}

### set group.multiplier
# sqrt by length
k2 <- k^0.5

# 1/2 interactions
k2[5:140] <- k2[5:140]*0.5
k2
#################################################
# fit
#fit <- grpsurv(as.matrix(X), y, group, penalty = "grLasso", lambda = seq(0.000001, 0.1, length.out=100))

# check the coef in lambda = 0.05
#k <- coef(fit, lambda = 0.05)
#k[k!=0]

# lambda plot
#plot(fit)

# cv.fit for finding lambda
cvfit <- cv.grpsurv(as.matrix(X), y, as.integer(factor(group))-1,
                    penalty="grLasso",
                    lambda = rev(seq(0.00001, 0.1, length.out=100)),
                    cv=5,
                    seed = 58465141)

cvfit$lambda
cvfit$lambda.min
lambda_min <- which(cvfit$lambda == cvfit$lambda.min)
# plot for cvfit
plot(cvfit, xlim = c(-2, -7), ylim = c(10,18))

# check the coef in lambda = cvfit
kk <- coef(cvfit, lambda = cvfit$lambda.min)
kk <- kk[kk!=0]

sel <- unique(substr(names(cvfit$fit$beta[,lambda_min][cvfit$fit$beta[,lambda_min] != 0]),1,9))

fit <- grpsurv(as.matrix(X), y, as.integer(factor(group))-1,
               penalty = "grLasso",
               group.multiplier = k)

plot(fit)

#######################################################################################################################################
### Visualization

# selected example brain roi
selected_interaction <- names(b)[c(29,55)]

brainregion <- data.frame(region = selected_interaction,
           gene = as.factor("GBA"),
           stringsAsFactors = F)

brainregion %>%
  ggseg(atlas = dk,
        colour = "grey37",
        mapping = aes(fill = gene))

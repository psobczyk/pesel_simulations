#' 
#' Piotr.Sobczyk (at) pwr.edu.pl 
#' 
#' 
#' 
path_to_matlab="/usr/local/MATLAB/R2013a/bin/matlab"
 
source("data_generation_fixed_effects.R")

library(varclust)
library(FactoMineR)
library(pryr)
library(softImpute)

args <- commandArgs(trailingOnly = TRUE)
print(args)
if(identical(args, character(0))) args <- c("")
id=args[[1]]

SNRs <- c(0.25, 0.375, 0.5, 0.625, 0.75, 1, 2, 4, 8)
vars <- c(50, 75, 100, 150, 200, 250, 400, 800, 1600)
numb.repetitions <- 100
k <- 5
n <- 100
scale <- TRUE
maxPC <- 2*k

methods.names <- c("PESEL n hetero", "Laplace evidence", "PESEL p hetero",
									 "GCV", "PESEL p homo", "Passemier", "CSV")

results <- NULL

tabela <- matrix(0, nrow=length(methods.names), ncol=numb.repetitions)

for(i in 1:length(SNRs)){
  for(j in 1:length(vars)){
    SNR <- SNRs[i]
    var <- vars[j]
    sim.data_old <- data.simulation.model(n, SNR, var, k)
    for(rep in 1:numb.repetitions) {
      sim.data <- data.fixed.model(object=sim.data_old, n = n, SNR = SNR, 
                                   numb.vars = var, k = k, scale = scale)
      
      BICs <- sapply(1:maxPC, function(j) varclust:::pca.BIC(sim.data$X, j))
      tabela[1,rep] <- which.max(BICs)
      BICs <- sapply(1:maxPC, function(j) varclust:::pca.Laplace(sim.data$X, j))
      tabela[2,rep] <- which.max(BICs)
      BICs <- sapply(1:maxPC, function(j) varclust:::pca.new.BIC(sim.data$X, j))
      tabela[3,rep] <- which.max(BICs)
      estimatedNpc <- estim_ncp(sim.data$X, ncp.min = 0, ncp.max = maxPC, method =  "GCV")$ncp
      tabela[4,rep] <- estimatedNpc
      BICs <- sapply(1:maxPC, function(j) varclust:::rajan.BIC(sim.data$X, j))
      tabela[5,rep] <- which.max(BICs)
      
      xFileName = paste0("x_", id, ".csv")
      write.table(sim.data$X, sep=",", dec=".", file = xFileName, col.names=FALSE, row.names = FALSE)
      command = paste0(path_to_matlab, " -nodisplay -r \"id='", id, 
                       "'; kmax=", maxPC, "; run('estimating_pcs_passemier2.m'); exit\"")
      system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)
      resFileName = paste0("res_", id, ".csv")
      result = read.csv(resFileName)
      
      if(is.na(result$sigma)){
        tabela[6,rep] <- 1
      } else{
        tabela[6,rep] <- min(result$k,maxPC)
      }
      
      n.folds <- 10
      folds <- NULL
      p <- var
      fold.size <- n*p/n.folds
      free.indices <- 1:(n*p)
      selected <- NULL
      temp <- sample(free.indices, size = fold.size, replace = FALSE)
      folds[temp] <- 1
      selected <- c(selected, temp)
      for(fold.i in 2:n.folds){
        temp <- sample(free.indices[-selected], size = fold.size, replace = FALSE)
        folds[temp] <- fold.i
        selected <- c(selected, temp)
      }
      
      ### creating lambda grid
      lam0 <- lambda0(sim.data$X)
      lamseq=exp(seq(from=log(lam0),to=log(1),length=20))
      
      ### cross-validation
      error <- matrix(0, nrow=n.folds, ncol=20)
      for(fold.i in 1:n.folds){
        xs <- sim.data$X
        xs[which(folds==fold.i)] <- NA
        for(lambda.j in seq_along(lamseq)){
          lambda=lamseq[lambda.j]
          fits=softImpute(xs, rank.max=5, type="svd", lambda = lambda)
          if(fits$d[1]!=0){
            error[fold.i,lambda.j] <- sum((sim.data$X-fits$u%*%diag(fits$d)%*%t(fits$v))[which(folds==fold.i)]^2)
          } else {
            error[fold.i,lambda.j] <- sum((sim.data$X[which(folds==fold.i)])^2)
          }
        }
      }
      
      #estimating sigma
      fits=softImpute(sim.data$X, rank.max = 5, lambda = lamseq[which.min(colMeans(error))],type="svd")
      deg.free <- sum(fits$d>0)
      if(deg.free!=0){
        sigma.hat <- sum((sim.data$X-fits$u%*%diag(fits$d)%*%t(fits$v))^2)/(n*(p-deg.free))
      } else{
        sigma.hat <- sum((sim.data$X)^2)/(n*(p-deg.free))
      }
      
      command = paste0(path_to_matlab, " -nodisplay -r \"id='", id,
                       "'; maxPC=", maxPC, "; sigma=", sigma.hat, "; run('csv_function'); exit\"")
      system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)
      
      result = read.csv(resFileName)
      
      if(is.na(result$sigma)){
      	tabela[7,rep] <- 1
      } else{
      	tabela[7,rep] <- min(result$k,maxPC)
      }
      results[[paste0(SNR, "_", var)]] <- tabela
      
    }
  }
}

# filename <- paste0("fixed_effects_simulations_", args[[1]], "_", gsub("-", "_", Sys.Date()), ".Rdata")
# save.image(filename)

####################
####### Plot #######
####################


library(dplyr)
library(ggplot2)
background="white"
main_color="darkred"
PeselPlotTheme <- theme(
	legend.position = "bottom",
	legend.direction = "horizontal",
	legend.box = "horizontal",
	legend.title = element_text(face = "italic", size = 17),
	legend.background = element_rect(fill = background),
	legend.key = element_rect(fill = background, colour = background),
	legend.text = element_text(size = 16),
	plot.background = element_rect(fill = background, colour = background),
	panel.background = element_rect(fill = background),
	panel.background = element_rect(fill = "white"),
	plot.title = element_text(face = "italic", size = 24, vjust = 1),
	panel.grid.major.x = element_blank(),
	panel.grid.major.y = element_line(),
	panel.grid.minor.y = element_line(),
	panel.grid.minor.x = element_line(colour = "gray95", size=1.5),
	axis.title = element_text(size = 17),
	axis.text = element_text(size = 15),
	strip.text.x = element_text(size = 17, colour = "white"),
	strip.background = element_rect(fill = main_color),
	axis.ticks.x = element_blank()
)

y <- NULL
result <- results
characteristics <- lapply(strsplit(names(result), "_"), as.numeric)
for(i in 1:length(result)){
	y <- rbind(y, data.frame(pc=as.vector(result[[i]]), method=rep(method.names, numb.repetitions),
													 snr=characteristics[[i]][1], variables=characteristics[[i]][2]))
}
y$pc <- pmax(1,floor(y$pc))


selected.snrs <- c(0.250, 0.375, 0.500, 0.625, 0.750, 1.000, 2, 4, 8)
selected.variables.number <- c(150, 1600)

selected.method2 <- method.names[c(1,2,3,4,6,7)]

y %>% filter(method %in% selected.method2,
						 snr %in% selected.snrs,
						 variables %in% selected.variables.number) %>% 
	mutate(snr=factor(snr)) %>% group_by(snr, variables, method) %>% 
	summarise(MeanPCs=mean(pc), SdPCs=sd(pc)/sqrt(n())) -> data5

new.lables=expression(paste("PESEL", phantom()["n"]^{"hetero"}), "CSV", "GCV", 
											paste("PESEL", phantom()["p"]^{"hetero"}), "Laplace Evidence", 
											"Passemier")

ggplot(data5, aes(x=snr, y=MeanPCs)) + 
	geom_line(aes(group=method, color=method, linetype=method), size=1.2) +
	# geom_errorbar(aes(ymin=MeanPCs-SdPCs, ymax=MeanPCs+SdPCs, color=method, alpha=1)) +
	facet_grid(.~variables) +
	ggtitle(paste0("Estimated number of PCs as a function of SNR")) +
	xlab("Signal to noise ratio (SNR)") +
	ylab("Mean estimated number of principal components") +
	scale_linetype_manual(values=c("solid", "dashed", "dotdash", "dotted", "longdash", "twodash")[seq_along(method.names)],
												labels = new.lables) +
	scale_color_discrete("Method", labels = new.lables) +
	scale_x_discrete(breaks=selected.snrs[-c(2,4)]) +
	PeselPlotTheme +
	guides(colour=guide_legend(nrow=2,byrow=TRUE), alpha=FALSE,
				 linetype=guide_legend("Method",nrow=2,byrow=TRUE, override.aes = c(size=1.7), 
				 											label.position = "bottom", label.hjust = 0.5,
				 											title.position = "left", keyheight=3, keywidth=14))

ggsave(filename = "figure3.png", height = 8, width = 14)
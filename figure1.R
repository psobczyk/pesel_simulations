
library(ggplot2)
# library(devtools)
# install_github("psobczyk/varclust")
library(varclust)

simulation_homogenous_eigenvalues <- function (SNR, k, n, numb.vars, scale = TRUE) {
  sigma <- 1/SNR
  X <- NULL 
  SIGNAL <- replicate(numb.vars, rnorm(n, 0, 1))
  
  for(i in 1:7){
  	pc.decomp <- svd(SIGNAL)
  	pc.decomp$d <- rep(mean(pc.decomp$d[1:k]), length(pc.decomp$d))
  	SIGNAL <- pc.decomp$u[,1:k] %*%  diag(pc.decomp$d[1:k]) %*%  t(pc.decomp$v[,1:k])
  	SIGNAL <- scale(SIGNAL)
  }
  X <- SIGNAL + replicate(numb.vars, sigma*rnorm(n))
  if(scale) X = scale(X)
  X
}

run_homogenous_simulation_for_given_sigma <- function (SNR, k, n, numb.vars, runMax = 10) {
  homoMean <- 0
  notHomoMean <- 0
  for(i in 1:runMax){
    X <- simulation_homogenous_eigenvalues(SNR, k, n, numb.vars, scale = TRUE)
    homoMean <- homoMean + pesel::pesel(X, method = "homo")
    notHomoMean <- notHomoMean + pesel::pesel(X, method = "hetero")
  }
  c(homoMean/runMax, notHomoMean/runMax)
}

simulation_non_homogenous_eigenvalues <- function (SNR, k, n, numb.vars, scale = TRUE) {
  sigma <- 1/SNR
  X <- NULL 

  SIGNAL <- replicate(numb.vars, rnorm(n, 0, 1))
  pc.decomp <- svd(SIGNAL)
  pc.decomp$d <- sum(pc.decomp$d[1:k])*2^(-1:(-length(pc.decomp$d)))
  SIGNAL <- pc.decomp$u[,1:k] %*%  diag(pc.decomp$d[1:k]) %*%  t(pc.decomp$v[,1:k])
  
  SIGNAL <- scale(SIGNAL)
  X <- SIGNAL + replicate(numb.vars, sigma*rnorm(n))
  if(scale){
    X = scale(X)
  }
  X
}

run_non_homogenous_simulation_for_given_sigma <- function (SNR, k, n, numb.vars, runMax = 10) {
  homoMean <- 0
  notHomoMean <- 0
  for(i in 1:runMax){
    X <- simulation_non_homogenous_eigenvalues(SNR, k, n, numb.vars, scale = TRUE)
    homoMean <- homoMean + pesel::pesel(X, method = "homo")
    notHomoMean <- notHomoMean + pesel::pesel(X, method = "hetero")
  }
  c(homoMean/runMax, notHomoMean/runMax)
}


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


n = 100
numb.vars = 50
k = 5
numb.runs <- 1000
SNRs <- seq(0.6, 0.7, 0.005)
dane_symulacje <- sapply(SNRs, function(x) run_homogenous_simulation_for_given_sigma(x, k, n, numb.vars, numb.runs))
dane_symulacje <- data.frame(rbind(cbind(SNR=SNRs, meanPC=dane_symulacje[1,], method="homo"),
																	 cbind(SNR=SNRs, meanPC=dane_symulacje[2,], method="nonHomo")))

SNRs <- seq(0.5, 3.5, 0.1)
dane_symulacje2 <- sapply(SNRs, function(x) run_non_homogenous_simulation_for_given_sigma(x, k, n, numb.vars, numb.runs))
dane_symulacje2 <- data.frame(rbind(cbind(SNR=SNRs, meanPC=dane_symulacje2[1,], method="homo"),
																		cbind(SNR=SNRs, meanPC=dane_symulacje2[2,], method="nonHomo")))

dane_symulacje3 <- cbind(dane_symulacje, type="All singular values equal")
dane_symulacje3 <- rbind(dane_symulacje3, cbind(dane_symulacje2, type="Singular values not equal"))
dane_symulacje3$meanPC <- as.numeric(as.character(dane_symulacje3$meanPC))
dane_symulacje3$SNR <- as.numeric(as.character(dane_symulacje3$SNR))

new.lables=expression(paste("PESEL", phantom()["n"]^{"homo"}),
											paste("PESEL", phantom()["n"]^{"hetero"}))

ggplot(dane_symulacje3) + geom_line(aes(x=SNR, y=meanPC, linetype=method,
																				group=method, color=method), size=1.2) +
	ggtitle("Comparing criteria derived for two different priors") +
	ylab("Mean estimated number of Principal Components") +
	facet_grid(.~type, scales="free") +
	scale_y_continuous(breaks=seq(1, 5, 0.5)) +
	scale_linetype_manual("Criterion", values=c("solid", "dashed"),labels=new.lables) +
	scale_color_manual("Criterion", values = c("blue", "red"), labels=new.lables) +
	sobczykPlotTheme +
	guides(linetype=guide_legend(title = "Criterion", override.aes = c(size=1.4), 
															 label.position = "bottom", label.hjust = 0.5,
															 title.position = "left", keyheight=3, keywidth=12))

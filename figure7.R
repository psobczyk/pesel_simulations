#' 
#' Piotr.Sobczyk (at) pwr.edu.pl 
#' 
#' 

library(varclust)
library(FactoMineR)
library(xtable)
library(MetabolAnalyze)

data(UrineSpectra)

#BIC for PPCA, PESEL_n
mdlfit<-ppca.metabol(UrineSpectra[[1]], minq=1, maxq=3, scale="none", plot.BIC = TRUE, printout = FALSE)

#GCV
estim_ncp(UrineSpectra[[1]], ncp.min = 1, ncp.max = 10, scale = FALSE)$ncp


# PESEL_p
varclust::estim.npc(UrineSpectra[[1]])
#To fit we use subspace clustering method
mlcc.fit <- varclust::mlcc.bic((UrineSpectra[[1]]), numb.clusters = 1, numb.runs = 1, max.dim = 8)

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

data <- data.frame(id=1:18, score=mlcc.fit$factors[[1]], groups=as.factor(group))
ggplot(data, aes(x=id,y=X1.1, col=groups, shape=groups)) + geom_point(size=9) + 
  geom_abline(slope=0,intercept=-1.95) +
  xlab("Index") + ylab("Factor scores") +
  scale_color_manual("Groups", values=c("red", "blue")) +
  scale_shape_discrete("Groups") +
  PeselPlotTheme
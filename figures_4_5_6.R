#' 
#' Piotr.Sobczyk (at) pwr.edu.pl 
#' 
#' 
source("data_generation_fixed_effects.R")

library(pryr) #for partial function evaluation
library(varclust)
library(FactoMineR)
library(softImpute)

args <- commandArgs(trailingOnly = TRUE)

if(identical(args, character(0))) args <- c("")
id=args[[1]]
#set seed based on ID
set.seed(strtoi(charToRaw(id)))

SNRs <- c(0.25, 0.375, 0.5, 0.625, 0.75, 1, 4, 8)
vars <- c(50, 75, 150, 400, 800)
n <- 100
numb.repetitions <- 100

#probabilistic assumption violation parameters
degrees.freedom <- 3
additional.ratio <- 0.5

results <- list()

print("starting...")

for(k in c(5)){
  results[[paste0("student.noise_df_", degrees.freedom, "_", k)]] =
    compare_methods(partial(data.simulation.student.noise, df=degrees.freedom), numb.repetitions = numb.repetitions,
                    n = n, SNRs = SNRs, vars = vars, k = k, scale = TRUE, id = args[[1]])
  print("student.noise done")
  
  results[[paste0("additional.variables_ratio_", additional.ratio, "_", k)]] =
    compare_methods(partial(data.simulation.additional.variables, ratio = additional.ratio),
                                                     numb.repetitions = numb.repetitions, n = n, SNRs = SNRs, vars = vars,
                                                     k = k, scale = TRUE, id = args[[1]])
  print("additional.variables done")
}

# filename <- paste0("robustness_simulations_", args[[1]], "_", gsub("-", "_", Sys.Date()), ".Rdata")
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

scheme.names <- c("Student noise", "Additional variables")
methods.names <- c("PESEL n hetero", "Laplace evidence", "PESEL p hetero",
									 "GCV", "PESEL p homo", "Passemier", "CSV")


y <- NULL

for(sim.scenario in seq_along(scheme.names)){
	result <- results[[sim.scenario]]
	characteristics <- lapply(strsplit(names(result), "_"), as.numeric)
	for(i in 1:length(result)){
		y <- rbind(y, data.frame(pc=as.vector(result[[i]]), method=rep(method.names, numb.repetitions),
														 snr=characteristics[[i]][1], variables=characteristics[[i]][2],
														 scenario=scheme.names[sim.scenario]))
	}
}
y$pc <- pmax(1,floor(y$pc))


###############################################
################ Figure 4 #####################
###############################################

selected.variables.number <- c(150, 800)
selected.snrs <- c(0.375, 0.625, 0.75, 1, 2, 8)
sim.scenario = 1
selected.method <- method.names[c(1,3,5,7)]

y %>% filter(scenario==scheme.names[sim.scenario], method==selected.method[1]) %>% 
	group_by(snr, variables, method) %>%
	summarise(count=n()) %>% ungroup %>% select(count) %>% unlist %>% head(1) -> numb.repetitions
data4 <- y %>% filter(scenario==scheme.names[sim.scenario], 
											method %in% selected.method,
											variables %in% selected.variables.number,
											snr %in% selected.snrs) %>% 
	mutate(snr=factor(snr)) %>%
	group_by(snr, variables, method, pc) %>%
	summarise(count=n())
y %>% filter(scenario==scheme.names[sim.scenario], 
						 method %in% selected.method,
						 variables %in% selected.variables.number,
						 snr %in% selected.snrs) %>% 
	mutate(snr=factor(snr)) %>% group_by(snr, variables, method) %>% 
	summarise(MeanPCs=mean(pc)) -> data5


title <- paste0(scheme.names[sim.scenario], ". Estimated number of PCs as a function of SNR.")
new.lables=expression(paste("PESEL", phantom()["n"]^{"hetero"}), "GCV", 
											paste("PESEL", phantom()["p"]^{"hetero"}), "Passemier")

ggplot(data4, aes(x=snr, colour=method, group=method, shape=method)) +
	geom_line(data=data5, aes(x=snr, y=MeanPCs)) +
	geom_label(aes(y=pc-0.4, label=count/numb.repetitions, fill=method), color="white", position = position_dodge(width = 0.70)) +
	geom_point(aes(y=pc,size=count/numb.repetitions), position = position_dodge(width = 0.70)) +
	facet_grid(.~variables) +
	ylab("Number of PCs") +
	xlab("SNR") +
	scale_y_continuous(breaks=pretty(data4$pc)) +
	scale_size("Frequency", range = c(1,10)) +
	scale_color_discrete("Method", labels = new.lables) +
	scale_shape_manual("Method", values = 15:18, labels = new.lables) +
	guides(fill=FALSE, label=FALSE,shape=guide_legend(override.aes = list(size=10, linetype=0),
																										label.position = "bottom", 
																										title.position = "left",
																										title.hjust = 1)) +
	guides(size=guide_legend(label.position = "bottom", title.position = "left", keyheight=4, keywidth=3)) +
	geom_vline(xintercept=seq(1.5, length(unique(data4$snr))-0.5, 1), 
						 lwd=1, colour="gray80", linetype=2) +
	ggtitle(title) +
	PeselPlotTheme

ggsave(filename = "figure4.png", height = 8, width = 14)

###############################################
################ Figure 5 #####################
###############################################

selected.variables.number <- c(150, 800)
selected.snrs <- c(0.25, 0.5, 0.63, 0.75, 1, 2.00, 8)
sim.scenario = 2
selected.method <- method.names[c(1,3,5,7)]

y %>% filter(scenario==scheme.names[sim.scenario], method==selected.method[1]) %>% 
	group_by(snr, variables, method) %>%
	summarise(count=n()) %>% ungroup %>% select(count) %>% unlist %>% head(1) -> numb.repetitions
data4 <- y %>% filter(scenario==scheme.names[sim.scenario], 
											method %in% selected.method,
											variables %in% selected.variables.number,
											snr %in% selected.snrs) %>% 
	mutate(snr=factor(snr)) %>%
	group_by(snr, variables, method, pc) %>%
	summarise(count=n())
y %>% filter(scenario==scheme.names[sim.scenario], 
						 method %in% selected.method,
						 variables %in% selected.variables.number,
						 snr %in% selected.snrs) %>% 
	mutate(snr=factor(snr)) %>% group_by(snr, variables, method) %>% 
	summarise(MeanPCs=mean(pc)) -> data5


title <- paste0(scheme.names[sim.scenario], ". Estimated number of PCs as a function of SNR.")
new.lables=expression(paste("PESEL", phantom()["n"]^{"hetero"}), "GCV", 
											paste("PESEL", phantom()["p"]^{"hetero"}), "Passemier")

ggplot(data4, aes(x=snr, colour=method, group=method, shape=method)) +
	geom_line(data=data5, aes(x=snr, y=MeanPCs)) +
	geom_label(aes(y=pc-0.4, label=count/numb.repetitions, fill=method), color="white", position = position_dodge(width = 0.70)) +
	geom_point(aes(y=pc,size=count/numb.repetitions), position = position_dodge(width = 0.70)) +
	facet_grid(.~variables) +
	ylab("Number of PCs") +
	xlab("SNR") +
	scale_y_continuous(breaks=pretty(data4$pc)) +
	scale_size("Frequency", range = c(1,10)) +
	scale_color_discrete("Method", labels = new.lables) +
	scale_shape_manual("Method", values = 15:18, labels = new.lables) +
	guides(fill=FALSE, label=FALSE,shape=guide_legend(override.aes = list(size=10, linetype=0),
																										label.position = "bottom", 
																										title.position = "left",
																										title.hjust = 1)) +
	guides(size=guide_legend(label.position = "bottom", title.position = "left", keyheight=4, keywidth=3)) +
	geom_vline(xintercept=seq(1.5, length(unique(data4$snr))-0.5, 1), 
						 lwd=1, colour="gray80", linetype=2) +
	ggtitle(title) +
	PeselPlotTheme

ggsave(filename = "figure5.png", height = 8, width = 14)

###############################################
################ Figure 6 #####################
###############################################
y$variables <- as.factor(y$variables)

selected.variables.number <- c(50, 75, 150, 400, 800)
selected.snrs <- c(1,4)
sim.scenario = 1
selected.method <- method.names[c(1,3,5,7)]

y %>% filter(scenario==scheme.names[sim.scenario], method==selected.method[1]) %>% 
	group_by(snr, variables, method) %>%
	summarise(count=n()) %>% ungroup %>% select(count) %>% unlist %>% head(1) -> numb.repetitions
data4 <- y %>% filter(scenario==scheme.names[sim.scenario], 
											method %in% selected.method,
											variables %in% selected.variables.number,
											snr %in% selected.snrs) %>% 
	mutate(snr=factor(snr)) %>%
	group_by(snr, variables, method, pc) %>%
	summarise(count=n())
y %>% filter(scenario==scheme.names[sim.scenario], 
						 method %in% selected.method,
						 variables %in% selected.variables.number,
						 snr %in% selected.snrs) %>% 
	mutate(snr=factor(snr)) %>% group_by(snr, variables, method) %>% 
	summarise(MeanPCs=mean(pc)) -> data5


title <- paste0(scheme.names[sim.scenario], ". Estimated number of PCs as a function of number of variables.")
new.lables=expression(paste("PESEL", phantom()["n"]^{"hetero"}), "GCV", 
											paste("PESEL", phantom()["p"]^{"hetero"}), "Passemier")

ggplot(data4, aes(x=variables, colour=method, group=method, shape=method)) +
	geom_line(data=data5, aes(x=variables, y=MeanPCs)) +
	geom_label(aes(y=pc-0.4, label=count/numb.repetitions, fill=method), color="white", position = position_dodge(width = 0.70)) +
	geom_point(aes(y=pc,size=count/numb.repetitions), position = position_dodge(width = 0.70)) +
	facet_grid(.~snr) +
	ylab("Number of PCs") +
	xlab("Variables") +
	scale_y_continuous(breaks=pretty(data4$pc)) +
	scale_size("Frequency", range = c(1,10)) +
	scale_color_discrete("Method", labels = new.lables) +
	scale_shape_manual("Method", values = 15:18, labels = new.lables) +
	guides(fill=FALSE, label=FALSE,shape=guide_legend(override.aes = list(size=10, linetype=0),
																										label.position = "bottom", 
																										title.position = "left",
																										title.hjust = 1)) +
	guides(size=guide_legend(label.position = "bottom", title.position = "left", keyheight=4, keywidth=3)) +
	geom_vline(xintercept=seq(1.5, length(unique(data4$variables))-0.5, 1), 
						 lwd=1, colour="gray80", linetype=2) +
	ggtitle(title) +
	PeselPlotTheme

ggsave(filename = "figure6.png", height = 8, width = 14)
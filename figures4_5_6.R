
library(dplyr)

source("plot_templates.R")

date = NULL
if(is.null(date)) stop("Please specify variable date in file 'figures_4_5_6.R'")

scheme.names <- c("Student noise", "Additional variables", "Lognormal noise")
method.names <- c("Minka BIC", "Laplace evidence", "Properly dervied BIC",
                  "GCV", "Rajan - penalized version", "Passemier", "CSV")

y <- NULL
for(batch in c("A", "B", "C", "D", "E")){
  f <- paste0("data/robustness_fixed_simulations", batch, "_n_100_", date, ".Rdata")
  load(f)
  
  for(sim.scenario in c(1,2,3)){
    result <- results[[sim.scenario]]
    characteristics <- lapply(strsplit(names(result), "_"), as.numeric)
    for(i in 1:length(result)){
      y <- rbind(y, data.frame(pc=as.vector(result[[i]]), method=rep(method.names, numb.repetitions),
                               snr=characteristics[[i]][1], variables=characteristics[[i]][2],
                               scenario=scheme.names[sim.scenario]))
    }
  }
}
y$method <- factor(y$method, levels = method.names[c(1:2,4,3,5:7)])

###############################################
############## VARS PLOTS #####################
###############################################

selected.variables.number <- c(150, 800)
selected.snrs <- c(0.25, 0.5, 0.75, 1, 2, 4)
sim.scenario = 1
selected.method <- method.names[c(1,3,4,7)]

selected.variables.number <- c(150, 800)
selected.snrs <- c(0.25, 0.5, 0.63, 0.75, 1, 2.00, 8)
sim.scenario = 2
selected.method <- method.names[c(1,3,4,7)]

selected.variables.number <- c(150, 800)
selected.snrs <- c(0.25, 0.375, 0.5, 0.63, 1, 2.00, 8)
sim.scenario = 3
selected.method <- method.names[c(1,3,4,7)]


plot_facet_vars(y, scheme.names, sim.scenario, selected.method, selected.snrs, selected.variables.number)
plot_name <- paste0(gsub(" ", "_", scheme.names[sim.scenario]), "_", 
                    paste0(selected.variables.number, collapse = "_"), ".png")
ggsave(filename = plot_name, height = 8, width = 16)


###############################################
############### SNR PLOTS #####################
###############################################


selected.variables.number <- c(50, 75, 150, 400, 800)
selected.snrs <- c(1,4)
sim.scenario = 1
selected.method <- method.names[c(1,3,4,7)]

plot_facet_snr(y, scheme.names, sim.scenario, selected.method, selected.snrs, selected.variables.number)
plot_name <- paste0(gsub(" ", "_", scheme.names[sim.scenario]), "_", 
                    paste0(selected.snrs, collapse = "_"), ".png")
ggsave(filename = plot_name, height = 8, width = 16)

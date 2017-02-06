
library(ggplot2)
library(dplyr)

background="white"
main_color="darkred"
sobczykPlotTheme <- theme(
  legend.position = "bottom",
  legend.direction = "horizontal",
  legend.box = "horizontal",
  legend.title = element_text(face = "italic", size = 17),
  legend.background = element_rect(fill = background),
  legend.key = element_rect(fill = background, colour = background),
  legend.text = element_text(size = 16),
  plot.background = element_rect(fill = background, colour = background),
  panel.background = element_rect(fill = background),
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

#' Create a custom facet labeller function
facet_labeller_snr <- function(variable, value) {
  paste0("SNR=", value)
}
facet_labeller_vars <- function(variable, value) {
  paste0("n=100, ", "p=", value)
}


plot_facet_vars <- function(y, scheme.names, sim.scenario, selected.method, selected.snrs, selected.variables.number){
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
    facet_grid(.~variables, labeller = facet_labeller_vars) +
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
    sobczykPlotTheme
}


plot_facet_snr <- function(y, scheme.names, sim.scenario, selected.method, selected.snrs, selected.variables.number){
  
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
  
  data4$variables <- as.factor(data4$variables)
  data5$variables <- as.factor(data5$variables)
  
  title <- paste0(scheme.names[sim.scenario], ". Estimated number of PCs as a function of number of variables.")
  new.lables=expression(paste("PESEL", phantom()["n"]^{"hetero"}), "GCV", 
                        paste("PESEL", phantom()["p"]^{"hetero"}), "Passemier")

  
  ggplot(data4, aes(x=variables, colour=method, group=method, shape=method)) +
    geom_line(data=data5, aes(x=variables, y=MeanPCs)) +
    geom_label(aes(y=pc-0.4, label=count/numb.repetitions, fill=method), color="white", position = position_dodge(width = 0.70)) +
    geom_point(aes(y=pc,size=count/numb.repetitions), position = position_dodge(width = 0.70)) +
    facet_grid(.~snr, labeller = facet_labeller_snr) +
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
    sobczykPlotTheme
}

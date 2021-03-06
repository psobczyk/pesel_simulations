#' 
#' Piotr Sobczyk, PWr 
#' 
#' We run this simulations parallel with 5 different arguments passed to the script (A, B, C, D and E)
#' 

path_to_matlab=NULL
if(is.null(path_to_matlab)) stop("Please specify variable path_to_matlab in file 'data_generation_for_figures_4_5_6.R'")

source("data_generation_robustness_fixed_effects.R")

library(pryr) #for partial function evaluation
library(pesel)
library(FactoMineR)
library(denoiseR)
library(softImpute)

args <- commandArgs(trailingOnly = TRUE)
print(args)
if(identical(args, character(0))) args <- c("")
id=args[[1]]

#set seed based on ID
set.seed(strtoi(charToRaw(id)))

SNRs <- c(0.25, 0.375, 0.5, 0.625, 0.75, 1, 4, 8)
vars <- c(50, 75, 150, 400, 800)
n <- 100
numb.repetitions <- 20

#this setting is just for veryfing that program works!!!!
# SNRs <- 2^(0:1)
# vars <- c(40,50)
# n <- 30
# numb.repetitions <- 2


#violation parameters
degrees.freedom = 3
additional.ratio = 0.5
mulog = 2
sdlog = 1.2

results <- list()

print("starting...")

for(k in c(5)){
  results[[paste0("student.noise_df_", degrees.freedom, "_", k)]] =
    compare_methods(partial(data.simulation.student.noise, df=degrees.freedom), numb.repetitions = numb.repetitions,
                    n = n, SNRs = SNRs, vars = vars, k = k, scale = TRUE, id = args[[1]], pathToMatlab = path_to_matlab)
  print("student noise done")

  results[[paste0("additional.variables_ratio_", additional.ratio, "_", k)]] =
    compare_methods(partial(data.simulation.additional.variables, ratio = additional.ratio),
                                                     numb.repetitions = numb.repetitions, n = n, SNRs = SNRs, vars = vars,
                                                     k = k, scale = TRUE, id = args[[1]], pathToMatlab = path_to_matlab)
  print("additional.variables done")
  
  results[[paste0("lognormal.noise_mu_", mulog, "_sd_", sdlog, "_", k)]] =
    compare_methods(partial(data.simulation.lognormal.noise, mu = mulog, sd = sdlog), numb.repetitions = numb.repetitions,
                    n = n, SNRs = SNRs, vars = vars, k = k, scale = TRUE, id = args[[1]], pathToMatlab = path_to_matlab)
  print("log-normal noise done")
}

filename <- paste0("data/robustness_fixed_simulations", args[[1]], "_n_", n, "_", gsub("-", "_", Sys.Date()), ".Rdata")
save.image(filename)

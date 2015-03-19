#########################################################################################
# Main program to run calculate.weights, calculate.error and w.protein.ratio functions.
# Lina Hultin Rosenberg 20140903 
# lina.hultin-rosenberg@scilifelab.se
#########################################################################################

#!/usr/bin/env Rscript

#Load dependencies and functions
suppressMessages(library(Hmisc))

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

source(file.path(script.basename, 'norm_sample_median.r'))
source(file.path(script.basename, 'getWeight.r'))
source(file.path(script.basename, 'calculate_weights.r'))
source(file.path(script.basename, 'calculate_error.r'))
source(file.path(script.basename, 'w_protein_ratio.r'))

args = commandArgs(TRUE)

if (length(args)<9) {
  stop("Too few input arguments. Provide the following input arguments:
       
       The name of the dataset to analyse

       The name of the base output file to be written to the current working directory
       This script will write 2 files using this base name

       The minimum intensity to include (noise level, must be higher than 0)
       
       Index for column with protein group accessions
       
       Index for column with proteins in protein group
       
       Indices (a,b,c,etc) for columns with quantitative data
       
       Index for column to use as numerator in weight calculation
       
       Index for column to use as denominator in weight calculation
       
       Index for column to use as denominator in final protein ratio calculation
       
       Optional index for column to use as denominator in final protein ratio calculation (mean of two columns)
       
       ")
}

dataset = as.character(args[1])
outfilename = as.character(args[2])
quant.min = as.numeric(args[3])
group.index = as.numeric(args[4])
protein.index = as.numeric(args[5])
quant.index = as.integer(strsplit(args[6], ",")[[1]])
num.1 = as.numeric(args[7])
den.1 = as.numeric(args[8])
if (length(args)>9) {
  den.2 = c(as.numeric(args[9]),as.numeric(args[10]))
} else {
  den.2 = as.numeric(args[9])
}

#Call function to normalise peptides to same sample median
normalized_peptides = norm.sample.median(dataset, quant.index)

#Call function to calculate weights based on internal training set (technical duplicates)
weight_matrix <- calculate.weights(normalized_peptides, quant.min, num.1, den.1)

#Call function to calculate error on protein lelvel based on internal training set (technical duplicates)
weight_results <- calculate.error(normalized_peptides, weight_matrix, quant.min, group.index, num.1, den.1)

#Call function to calculate weighted protein quant for all quant columns
w.protein.ratio(normalized_peptides, outfilename, weight_matrix, weight_results,  quant.min, group.index, protein.index, quant.index, den.2)

#########################################################################################
# Main program to run calculate.weights, calculate.error and w.protein.ratio functions.
# Lina Hultin Rosenberg 20140903 
# lina.hultin-rosenberg@scilifelab.se
#########################################################################################

#!/usr/bin/env Rscript

#Load dependencies and functions
suppressMessages(library(Hmisc))

source('~/r_scripts/norm_sample_median.r')
source('~/r_scripts/getWeight.r')
source('~/r_scripts/calculate_weights.r')
source('~/r_scripts/calculate_error.r')
source('~/r_scripts/w_protein_ratio.r')

args = commandArgs(TRUE)

if (length(args)<9) {
  stop("Too few input arguments. Provide the following input arguments:
       
       The name of the dataset to analyse

       The name of the output file to be written to the current working directory

       The minimum intensity to include (noise level)
       
       Index for column with protein group accessions
       
       Index for column with proteins in protein group
       
       Range of indices (a:b) for columns with quantitative data
       
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
quant.range = strsplit(args[6],":")
quant.index = c(quant.range[[1]][1]:quant.range[[1]][2])
num.1 = as.numeric(args[7])
den.1 = as.numeric(args[8])
if (length(args)>9) {
  den.2 = c(as.numeric(args[9]),as.numeric(args[10]))
} else {
  den.2 = as.numeric(args[9])
}

#Call function to normalise peptides to same sample median
normfile = norm.sample.median(dataset, quant.index)

#Call function to calculate weights based on internal training set (technical duplicates)
weight_file <- calculate.weights(normfile, quant.min, num.1, den.1)

#Call function to calculate error on protein lelvel based on internal training set (technical duplicates)
weight_results_file <- calculate.error(normfile, weight_file, quant.min, group.index, num.1, den.1)

#Call function to calculate weighted protein quant for all quant columns
w.protein.ratio(normfile, weight_results_file,  quant.min, group.index, protein.index, quant.index, den.2)


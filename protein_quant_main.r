##############################################
# Main program to calculate protein ratios
# Lina Hultin Rosenberg 201409018 
# lina.hultin-rosenberg@scilifelab.se
##############################################

#!/usr/bin/env Rscript

#Load function needed
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
source(file.path(script.basename, 'norm_sample_median.r'))
source(file.path(script.basename, 'protein_ratio.r'))

args = commandArgs(TRUE)

if (length(args)<7) {
  stop("Too few input arguments. Provide the following input arguments:
       
       The full path to the tsv file to run protein quant on
       
       The name of the output file which will be written to the working directory
       
       The minimum intensity to include (noise level, must be higher than 0)
       
       Index for column with protein group accessions
       
       Index for column with proteins in protein group
       
       Indices (a,b,c,etc) for columns with quantitative data
       
       One or more indices for columns to use as denominators in protein ratio calculation (mean of columns is taken)
       ")
}

dataset = as.character(args[1])
outfilename = as.character(args[2])
quant.min = as.numeric(args[3])
group.index = as.numeric(args[4])
protein.index = as.numeric(args[5])
quant.index = as.integer(strsplit(args[6], ",")[[1]])
if (length(args)>7) {
  den = as.numeric(args[7:len(args)])
} else {
  den = as.numeric(args[7])
}
#Call function to normalise peptides to same sample median
normalized_peptides = norm.sample.median(dataset, quant.index)

#Call function to calculate protein ratios
protein.ratio(normalized_peptides, outfilename, quant.min, group.index, protein.index, quant.index, den)

################

norm.sample.median = function(filename, quant.index)
{
  
  #Load peptide data
  pep.data = read.delim(filename,header=TRUE,check.names=FALSE,row.names=NULL,sep="\t")
  
  #Extract quantitative data
  ratios = pep.data[,quant.index]
  
  norm_ratios = ratios
  sample.median = c()
  for (i in 1:ncol(ratios)) {
    sample.median[i] = median(ratios[,i],na.rm=TRUE)
  }
  #median.mean = mean(sample.median)
  for (i in 1:ncol(ratios)) {
    norm_ratios[, i] = (ratios[, i]/sample.median[i])
  }
  
  #Return normalised PSM table
  pep.norm = pep.data
  pep.norm[,quant.index] = norm_ratios
  
  return(pep.norm)
}

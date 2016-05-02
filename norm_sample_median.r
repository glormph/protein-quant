################
divide_column = function(numerator, denominator){
  divided = numerator / denominator
  return(divided)
  }


norm.sample.median = function(filename, quant.min, quant.index, denominator)
{
  
  #Load peptide data
  pep.data = read.delim(filename,header=TRUE,check.names=FALSE,row.names=NULL,sep="\t")
  channels = colnames(pep.data)[quant.index] 

  # Set NA to channel if below noise level  
  pep.quant = pep.data[, quant.index] 
  pep.quant[pep.quant < quant.min] = NA
  pep.data[, quant.index] = pep.quant

  # Define denominator column
  if (length(denominator)==1) {
    quant.den = pep.data[,denominator]
  } else {
    quant.den = apply(cbind(pep.data[,denominator]),1,mean,na.rm=TRUE)  
  }
   
  #Extract quantitative data and ratios
  ratios = apply(pep.data[, quant.index], 2, divide_column, denominator=quant.den)
  norm_ratios = ratios
  sample.median = c()
  print("Sample medians per channel:")
  for (i in 1:ncol(ratios)) {
    sample.median[i] = median(ratios[,i],na.rm=TRUE)
    report = sprintf("channel %s: %f", channels[i], sample.median[i])
    print(report)
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

########################################################################
# Function to calculate weights based on peptide intensity and variance
# Lina Hultin Rosenberg 20140904 
# lina.hultin-rosenberg@scilifelab.se
########################################################################

calculate.weights = function(pep.data, quant.min, num, den) {
  
  #Define parameters for weight calculation
  ratio = 1           #Expected ratio between technical duplicates
  bins = 8            #Number of bins to divide peptide in for weight calculation
  
  ##Calculate weight-matrix for technical duplicates
  
  #Extract quant columns for weight calculation
  pep.quant = pep.data[,c(num,den)]
  
  #Exchange everything smaller than quant.min by NA
  pep.quant[pep.quant<quant.min] = NA
  
  #Remove missing data (peptides with NA in at least one quant column)
  index.na = unique(c(which(is.na(pep.quant[,1])),which(is.na(pep.quant[,2]))))
  if (length(index.na)>0) {
    pep.quant = pep.quant[-index.na,]
  } 
  
  #Calculate variance as deviance from expected ratio (on log2 scale)
  estimated.ratio = log2(pep.quant[,1]/pep.quant[,2])
  expected.ratio = log2(ratio)
  variance = (estimated.ratio - expected.ratio)^2
  
  #Calculate reference intensity as minimum intensity over quant columns (log2)
  reference = log2(apply(pep.quant,1,min))
  
  #Create variance matrix and sort in ascending order of reference intensity 
  variance.matrix = cbind(reference, variance) 
  variance.matrix = variance.matrix[order(variance.matrix[,1]),]
  
  #Create weight matrix (call function getWeight)
  weight.matrix = getWeightMatrix(variance.matrix,bins)
  
  return(weight.matrix)
  ##Plot ratios versus minimum intensity
  #filename = paste(folder.weights,dataset,"_ratio_plot.tif",sep="")
  #tiff(file=filename)
  #filename = paste(folder.weights,dataset,"_ratio_plot.png",sep="")
  #png(file=filename)
  #par(mar=c(2.5,2.5,0.5,0.5),mgp=c(1.4,0.5,0))
  #plot(reference,estimated.ratio,main="",xlab="log2 of minimum intensity",ylab="log2 of ratio",pch=20)
  #abline(a=0,b=0)
  #abline(v=weight.matrix[-nrow(weight.matrix),3],col=grey(0.8))
  #dev.off()
  #
  ###Plot weights for bins
  ##filename = paste(folder.weights,dataset,"_weight_plot.tif",sep="")
  ##tiff(file=filename)
  #filename = paste(folder.weights,dataset,"_weight_plot.png",sep="")
  #png(file=filename)
  #par(mar=c(2.5,2.5,1,0.5),mgp=c(1.4,0.5,0))
  #barplot(weight.matrix[,4],main="",names.arg="",xlab="Bins",ylab="Weights")
  #dev.off()
  
}
  


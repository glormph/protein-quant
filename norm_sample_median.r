#Function for normalising to same sample median
#Lina Hultin-Rosenberg 20120222

norm.sample.median = function(dataset,folder,quant.index)
{
  
  #Define folder to save output to
  folder.data = paste(folder,"data/",sep="")
  
  #Load peptide data
  filename = paste(folder.data,dataset,".txt",sep="")
  pep.data = read.delim(filename,header=TRUE,check.names=FALSE,row.names=NULL,sep="\t")
  
  #Extract quantitative data
  int = pep.data[,quant.index]
  
  int.norm = int
  sample.median = c()
  for (i in 1:ncol(int)) {
    sample.median[i] = median(int[,i],na.rm=TRUE)
  }
  median.mean = mean(sample.median)
  for (i in 1:ncol(int)) {
    int.norm[,i] = (int[,i]/sample.median[i])*median.mean
  }
  
  #Save normalised peptide data
  pep.norm = pep.data
  pep.norm[,quant.index] = int.norm
  filename = paste(folder.data,dataset,"_norm.txt",sep="")
  write.table(pep.norm,file=filename,row.names=FALSE,sep="\t")
  
}
#######################################################################################
# Function to calculate weighted protein quantities and errors for technical duplicates
# Lina Hultin Rosenberg 20140903 
# lina.hultin-rosenberg@scilifelab.se
#######################################################################################


calculate.error = function(dataset,folder,quant.min,group.index,num,den) {
  
  #Define parameters for weight calculation
  ratio = 1           #Expected ratio between duplicates
  bins = 8            #Number of bins for weight calculation
  
  #Define folders to save output to
  folder.data = paste(folder,"data/",sep="")
  folder.weights =  paste(folder,"weights/",sep="")
  
  #Load peptide data
  filename = paste(folder.data,dataset,".txt",sep="")
  pep = read.delim(filename,header=TRUE,check.names=FALSE,row.names=NULL,sep="\t")
  
  #Remove peptides without protein group accession and peptides shared between several protein groups
  pep = pep[pep[,group.index]!="",]
  pep = pep[-grep(";",pep[,group.index]),]
  
  #Extract quant columns
  pep.quant = pep[,c(num,den)]
  
  #Exchange everything smaller than quant.min by NA
  pep.quant[pep.quant<quant.min] = NA
 
  #Remove peptides (rows) with NA in at least one quant column
  index.na = unique(c(which(is.na(pep.quant[,1])),which(is.na(pep.quant[,2]))))
  if (length(index.na)>0) {
    pep.quant = pep.quant[-index.na,]
    pep = pep[-index.na,]
  } 
  
  #Get unique proteins
  proteins = as.character(pep[,group.index])
  index.duplicated = which(duplicated(proteins))
  unique.proteins = proteins[-index.duplicated]
  no.proteins = length(unique.proteins)
  
  ##Calculate weighted protein quantities for technical duplicates ratio
  
  #Load weight matrix
  filename = paste(folder.weights,dataset,"_weight_matrix.txt",sep="")
  weight.matrix = read.delim(filename,header=TRUE,sep="\t")
  
  #Select quant columns to calculate protein quant for
  quant.num = pep.quant[,1]
  quant.den = pep.quant[,2]
  
  #Calculate weighted protein quant and error for each protein
  weight.results = data.frame(ratios=rep(0,no.proteins),protein.weights=rep(0,no.proteins),peptides=rep(0,no.proteins),errors=rep(0,no.proteins),row.names=unique.proteins)
  
  for (protein.i in 1:no.proteins) {
    
    quant.num.protein = quant.num[proteins==unique.proteins[protein.i]]
    quant.den.protein = quant.den[proteins==unique.proteins[protein.i]]
    
    #Get weights for peptide ratio, use minimum intensity as reference (log2)
    reference = log2(apply(cbind(quant.num.protein,quant.den.protein),1,min,na.rm=TRUE))
    
    #Call function getWeight
    ratio.weights = getWeight(reference,weight.matrix[,3],weight.matrix[,4])
    weight.factor = ratio.weights[,2]
    
    #Calculate peptide ratios
    pep.ratio = log2(quant.num.protein/quant.den.protein)
    
    #Calculate total protein weight as average over peptide weights
    weight.results$protein.weights[protein.i] = mean(weight.factor)
    
    #Use Hmisc function to calculate weighted average (protein quantity)
    protein.ratio.w = wtd.mean(pep.ratio,weight.factor) 
    
    #Move back from log space to regular ratio
    weight.results$ratios[protein.i] = 2^protein.ratio.w
    
    #Calculate relative error for weighted protein quant
    error = abs(ratio-weight.results$ratios[protein.i])
    weight.results$errors[protein.i] = (error/ratio)*100
    
    #Number of peptides for protein
    weight.results$peptides[protein.i] = length(pep.ratio)
  }
  
  #Check and remove NA ratios
  index.na = which(is.na(weight.results$ratios))
  weight.results[index.na,] = NA
  
  #Save weight results
  filename = paste(folder.weights,dataset,"_weight_results.txt",sep="")
  write.table(weight.results,file=filename,col.names=NA,sep="\t")
  
}



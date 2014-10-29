########################################
# Function to calculate protein ratios
# Lina Hultin Rosenberg 20140918 
# lina.hultin-rosenberg@scilifelab.se
########################################


protein.ratio = function(dataset,folder,quant.min,group.index,protein.index,quant.index,den) {
  
  #Define parameters for weight calculation
  ratio = 1           #Expected ratio between duplicates
  bins = 8            #Number of bins for weight calculation
  
  #Define folders to save output to
  folder.data = paste(folder,"data/",sep="")
  folder.weights =  paste(folder,"weights/",sep="")
  folder.results = paste(folder,"quant/",sep="")
  
  #Load peptide data
  filename = paste(folder.data,dataset,".txt",sep="")
  pep = read.delim(filename,header=TRUE,check.names=FALSE,row.names=NULL,sep="\t")
  
  #Remove peptides without protein group accession and peptides shared between several protein groups
  pep = pep[pep[,group.index]!="",]
  pep = pep[-grep(";",pep[,group.index]),]
  
  #Extract quant columns
  pep.quant = pep[,quant.index]
  
  #Exchange everything smaller than quant.min by NA
  pep.quant[pep.quant<quant.min] = NA
  
  #Remove peptides (rows) with NA in all quant columns
  remove.index = c()
  for (i in 1:nrow(pep.quant)) {
    index.na = which(is.na(pep.quant[i,]))
    
    if (length(index.na)==ncol(pep.quant)) {
      remove.index = c(remove.index,i)
    }
  }
  if (length(remove.index)>0) {
    pep.quant = pep.quant[-remove.index,]
    pep = pep[-remove.index,]
  }
  
  #Define denominator for protein ratio calculation
  if (length(den)==1) {
    quant.den = pep[,den]
  } else {
    quant.den = apply(cbind(pep[,den[1]],pep[,den[2]]),1,mean,na.rm=TRUE)  
  }
  
  #Get unique proteins
  proteins = as.character(pep[,group.index])
  proteins.all = as.character(pep[,protein.index])
  index.duplicated = which(duplicated(proteins))
  unique.proteins = proteins[-index.duplicated]
  unique.proteins.all = proteins.all[-index.duplicated]
  no.proteins = length(unique.proteins)
  
  ##Calculate protein ratios (all columns)
  
  protein.ratios = matrix(data=NA,nrow=no.proteins,ncol=ncol(pep.quant))
  protein.peptides = matrix(data=NA,nrow=no.proteins,ncol=ncol(pep.quant))
  
  ratio.list = c()
  ratio.i = 1
  for (i in 1:ncol(pep.quant)) {
    
    #Select quant columns to calculate protein ratio for
    quant.num = pep.quant[,i]   
    
    #Get peptides belonging to each protein
    for (protein.i in 1:no.proteins) {
      
      quant.num.protein = quant.num[proteins==unique.proteins[protein.i]]
      quant.den.protein = quant.den[proteins==unique.proteins[protein.i]]
      
      #Calculate peptide ratios and remove empty
      pep.ratio = quant.num.protein/quant.den.protein
      index.na = which(is.na(pep.ratio))
      if (length(index.na)>0) {
        pep.ratio = pep.ratio[-index.na]
      } 
      
      #Calculate protein ratio as median over peptide ratios
      protein.ratios[protein.i,ratio.i] = median(pep.ratio,na.rm=TRUE)
      
      #Number of peptides used for protein
      protein.peptides[protein.i,ratio.i] = length(pep.ratio)
    }
    ratio.i = ratio.i+1
  }
  rownames(protein.ratios) = unique.proteins
  rownames(protein.peptides) = unique.proteins
  
  #Save table
  quant.table = cbind(protein.ratios,protein.peptides)
  quant.columns = c(1:ncol(pep.quant))
  column.names = c(quant.columns,paste(quant.columns,"no.peptides",sep=".")) 
  colnames(quant.table) = column.names
  
  filename = paste(folder.results,dataset,"_protein_ratios.txt",sep="")
  write.table(quant.table,file=filename,col.names=NA,sep="\t")
  
}


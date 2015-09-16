###############################################################################
# Function to calculate weighted protein quantities and errors for all channels
# Lina Hultin Rosenberg 20140903 
# lina.hultin-rosenberg@scilifelab.se
###############################################################################


w.protein.ratio = function(pep, filename, weight.matrix, weight.results, quant.min,group.index, protein.index, quant.index, den) {
  
  #Define parameters for weight calculation
  ratio = 1           #Expected ratio between duplicates
  bins = 8            #Number of bins for weight calculation
  
  #Define number of peptides for plotting
  no.peptides = c(1,2,3,4,5,6,10)
  no.peptides.text = c("1","2","3","4","5","6-10",">10")
  
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
    quant.den = apply(cbind(pep[,den]),1,mean,na.rm=TRUE)  
  }
  
  #Get unique proteins
  proteins = as.character(pep[,group.index])
  proteins.all = as.character(pep[,protein.index])
  index.duplicated = which(duplicated(proteins))
  unique.proteins = proteins[-index.duplicated]
  unique.proteins.all = proteins.all[-index.duplicated]
  no.proteins = length(unique.proteins)
  
  ##Calculate weighted protein ratios (all columns)
  
  protein.weights = matrix(data=NA,nrow=no.proteins,ncol=ncol(pep.quant))
  protein.ratios = matrix(data=NA,nrow=no.proteins,ncol=ncol(pep.quant))
  protein.peptides = matrix(data=NA,nrow=no.proteins,ncol=ncol(pep.quant))
  
  ratio.list = c()
  ratio.i = 1
  for (i in 1:ncol(pep.quant)) {
    
    #Extract quant column to calculate protein quant for
    quant.num = pep.quant[,i]   
    
    #Get peptides belonging to each protein
    for (protein.i in 1:no.proteins) {
      
      quant.num.protein = quant.num[proteins==unique.proteins[protein.i]]
      quant.den.protein = quant.den[proteins==unique.proteins[protein.i]]
      
      #Get weights for peptide ratio, use minimum intensity as reference (log2)
      reference = log2(apply(cbind(quant.num.protein,quant.den.protein),1,min,na.rm=TRUE))
      
      #Call function getWeight
      ratio.weights = getWeight(reference,weight.matrix[,3],weight.matrix[,4])
      weight.factor = ratio.weights[,2]
      
      #Calculate peptide ratios and remove empty
      pep.ratio = log2(quant.num.protein/quant.den.protein)
      index.na = which(is.na(pep.ratio))
      if (length(index.na)>0) {
        pep.ratio = pep.ratio[-index.na]
        weight.factor = weight.factor[-index.na]
      } 
      
      #Calculate total protein weight as average over peptide weights
      protein.weights[protein.i,ratio.i] = mean(weight.factor)
      
      #Use Hmisc function to calculate weighted average (protein quantity)
      protein.ratio.w = wtd.mean(pep.ratio,weight.factor) 
      
      #Move back from log space to regular ratio
      protein.ratios[protein.i,ratio.i] = 2^protein.ratio.w
      
      #Number of peptides for protein
      protein.peptides[protein.i,ratio.i] = length(pep.ratio)
    }
    ratio.i = ratio.i+1
  }
  rownames(protein.weights) = unique.proteins
  rownames(protein.ratios) = unique.proteins
  rownames(protein.peptides) = unique.proteins
  
  #Load data for training set
  
  #Get loess function for error versus average weight of proteins (95%) - for technical duplicates
  
  #Plot error versus average weight of proteins (95% loess smoothers)
  #filename = paste(folder.results,dataset,"_protein_error_weights.tif",sep="")
  #tiff(file=filename)
  #filename = paste(folder.results,dataset,"_protein_error_weights.png",sep="")
  #png(file=filename)
  #par(mar=c(3,3,0.5,0.5),mgp=c(1.5,0.5,0))
  
  #For each number of peptides
  protein.weights.ref = c()
  protein.errors.ref = c()
  loess.list = list()
  
  for (i in 1:length(no.peptides)) {
    if (i<7) {
      protein.list = rownames(weight.results[weight.results$peptides>=no.peptides[i]&weight.results$peptides<no.peptides[i+1],])
    } else {
      protein.list = rownames(weight.results[weight.results$peptides>no.peptides[i],])
    }
    protein.weights.ref = weight.results[protein.list,]$protein.weights
    protein.errors.ref = weight.results[protein.list,]$errors
    
    #Order the errors in order of the weights
    weights.order = order(protein.weights.ref)
    protein.weights.ref = protein.weights.ref[weights.order]
    protein.errors.ref = protein.errors.ref[weights.order]
    
    #Calculate median error
    #protein.errors.95 = c()
    protein.errors.median = c()
    if (i==1) {
      weights.unique = unique(protein.weights.ref)
      for (weight.unique in weights.unique) {
        protein.errors.sel = protein.errors.ref[which(protein.weights.ref==weight.unique)]
        #protein.errors.95 = c(protein.errors.95,quantile(protein.errors.sel,probs=seq(0,1,0.05),na.rm=TRUE)[20])
        protein.errors.median = c(protein.errors.median,median(protein.errors.sel,na.rm=TRUE))
      }
      protein.errors.runmed = protein.errors.median
      weight.breaks = weights.unique
    } else if (i==2)  {
      weights.unique = unique(protein.weights.ref)
      for (weight.unique in weights.unique) {
        protein.errors.sel = protein.errors.ref[which(protein.weights.ref==weight.unique)]
        #protein.errors.95 = c(protein.errors.95,quantile(protein.errors.sel,probs=seq(0,1,0.05),na.rm=TRUE)[20])
        protein.errors.median = c(protein.errors.median,median(protein.errors.sel,na.rm=TRUE))
      }
      protein.errors.runmed = runmed(protein.errors.median,11)
      weight.breaks = weights.unique
    } else {
      weight.breaks = quantile(protein.weights.ref,probs=seq(0,1,0.02))
      for (j in 1:(length(weight.breaks)-1)) {
        if (j<length(weight.breaks)-1) {
          error.interval = protein.errors.ref[protein.weights.ref>=weight.breaks[j]&protein.weights.ref<weight.breaks[j+1]]
        } else {
          error.interval = protein.errors.ref[protein.weights.ref>=weight.breaks[j]&protein.weights.ref<=weight.breaks[j+1]]
        }
        #protein.errors.95 = c(protein.errors.95,quantile(error.interval,probs=seq(0,1,0.05),na.rm=TRUE)[20])
        protein.errors.median = c(protein.errors.median,median(error.interval,na.rm=TRUE))
      }
      
      #Fill upp empty values
      index.na = which(is.na(protein.errors.median))
      
      if (length(index.na)>(length(protein.errors.median)/2)) {
        print(paste("Warning, too few proteins with",no.peptides[i],"peptides, error calculation not performed"))
        loess.list[i] = NA
        next
      }
      for (index in index.na) {
        error.median = NA
        k = 1
        while (is.na(error.median)) {
          index.start = max(c(1,(index-k)))
          index.stop = min(c(length(protein.errors.median),(index+k)))
          error.median = median(protein.errors.median[c(index.start:index.stop)],na.rm=TRUE)
          k = k+1
        }
        protein.errors.median[index] = error.median
      }
      protein.errors.runmed = runmed(protein.errors.median,11)
      weight.breaks = weight.breaks[-length(weight.breaks)]
    }
    
    #Fit a loess curve to error
    data.errors = data.frame(x=weight.breaks,y=protein.errors.runmed)
    lo.errors = loess(y~x,data.errors,span=0.5,degree=1)
    loess.list[i] = list(lo.errors)
    
  #  if (i==1) {
  #    plot(data.errors$x,predict(lo.errors,data.errors$x),xlab="Protein weight",ylab="Relative error (%)",type="l",col=i,ylim=c(0,max(protein.errors.runmed)))  
  #  } else {
  #    lines(data.errors$x,predict(lo.errors,data.errors$x),col=i)
    
  }
  #legend("topright",no.peptides.text,col=c(1:7),lty=1)
  #dev.off()
  
  #Create matrix to save relative error estimated by loess
  error.matrix = matrix(NA,nrow=nrow(protein.ratios),ncol=ncol(protein.ratios),byrow=FALSE,dimnames=list(rownames(protein.ratios)))
  
  #Get relative error for all ratios
  for (i in 1:nrow(protein.weights)) {
    for (j in 1:ncol(protein.weights)) {
      weight = as.numeric(protein.weights[i,j])
      peptides = as.numeric(protein.peptides[i,j])
      if (peptides<6) {
        peptides.index = peptides
      } else if (peptides>5 & peptides<11) {
        peptides.index = 6
      } else {
        peptides.index = 7
      }
      if (peptides.index==0 | weight==0) {
        error.matrix[i,j] = NA
      } else {
        if (length(loess.list[[peptides.index]])>1) {
          error.matrix[i,j] = predict(loess.list[[peptides.index]],weight) 
        }
      }
    }
  }
  
  #Round error and weight
  error.matrix = round(error.matrix)
  protein.weights = round(protein.weights)
  
  #Save table
  quant.table = cbind(rownames(protein.ratios), unique.proteins.all, protein.ratios, protein.peptides, protein.weights, error.matrix)
  quant.columns = colnames(pep.quant)
  column.names = c("Protein accession", "Proteins in group", quant.columns, paste(quant.columns, "# peptides", sep=" "), paste(quant.columns, "weight", sep=" "), paste(quant.columns, "rel. error", sep=" ")) 
  
  write.table(quant.table, file=filename, col.names=column.names, row.names=F, sep="\t", quote=F)
}
  
  

###############
###############


protein.ratio = function(pep, filename_out, quant.min, group.index, protein.index, quant.index, den) {
  
  #Remove peptides without protein group accession and peptides shared between several protein groups
  pep = pep[pep[,group.index]!="",]
  pep = pep[-grep(";",pep[,group.index]),]
  
  #Extract quant columns
  pep.quant = pep[,quant.index]
  
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
  # Sort peptide table by protein
  pep = pep[order(pep[,group.index]),]
  
  # Get amount of proteins left, create empty table
  unique_proteins = unique(pep[,group.index])
  protein_ratios = matrix(data=NA, nrow=length(unique_proteins), ncol=length(quant.index))
  protein_psmno = matrix(data=NA, nrow=length(unique_proteins), ncol=length(quant.index))
  
  # Loop sorted PSM table to get PSMs for each protein.
  protein_index = 1
  collected_ratios = matrix(ncol=ncol(protein_ratios))
  lastprotein = pep[1, group.index]
  for (i in 1:dim(pep)[1]){
    if (pep[i, group.index] != lastprotein){
       rownames(protein_ratios)[protein_index] = pep[i, group.index]
       for (chan in 1:ncol(collected_ratios)){
         protein_ratios[protein_index, chan] = median(collected_ratios[,chan])
         protein_psmno[protein_index, chan] = length(collected_ratios[complete.cases(collected_ratios[, chan]), chan])
         }
       protein_index =+ 1
       collected_ratios = matrix(ncol=ncol(protein_ratios))
       lastprotein = pep[i, group.index]
       }
    collected_ratios = rbind(collected_ratios, pep[i, quant.index]) 

  quant.table = cbind(rownames(protein_ratios), protein_ratios, protein_psmno)
  quant.columns = colnames(pep[, quant.index])
  column.names = c("Protein accession", quant.columns, paste(quant.columns,"# quanted PSMs",sep=" - "))
  
  write.table(quant.table, file=filename_out, row.names=F, col.names=column.names, sep="\t", quote=F)
}

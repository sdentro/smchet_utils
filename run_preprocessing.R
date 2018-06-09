#### this below should go into a new dpclust_smchet_utils repo and dpc pipeline could be updated



############################################################################################
# Preprocessing helpers as dpclust3p cannot directly work with the supplied VCF files for now
############################################################################################
#' Create a file with the SNV loci
createLociFile = function(vcfdat, outfile, chrom_col, pos_col, ref_col, alt_col) {
  loci = vcfdat[, c(chrom_col, pos_col, ref_col, alt_col)]
  loci_file = "loci.txt"
  write.table(loci, file=outfile, sep="\t", quote=F, row.names=F, col.names=F)
}

#' Convenience function that transforms counts for alt and ref alleles into alleleCounter output
mutwt2allelecounts = function(counts.alt, counts.ref, allele.alt, allele.ref) {
  output = array(0, c(length(allele.ref), 4))
  nucleotides = c("A", "C", "G", "T")
  # Propagate the alt allele counts
  nucleo.index = match(allele.alt, nucleotides)
  for (i in 1:nrow(output)) {
    output[i,nucleo.index[i]] = counts.alt[i]
  }
  
  # Propagate the reference allele counts
  nucleo.index = match(allele.ref, nucleotides)
  for (i in 1:nrow(output)) {
    output[i,nucleo.index[i]] = counts.ref[i]
  }
  return(output)
}

#' Dump a file with allele counts in alleleCounter output format
createAlleleCountsFile = function(vcfdat, datacol, namecol, outfile) {
  # Allelic depths for the ref and alt alleles (AD), Approximate read depth (DP)
  tumour_stat = data.frame(do.call(rbind, strsplit(as.vector(vcfdat[,datacol]), split = ":", fixed = TRUE)))
  colnames(tumour_stat) = strsplit(as.vector(unique(vcfdat[,namecol])),':')[[1]]
  # get the number of mutant reads into mutReads and the total number of reads at each mutation into totalReads, then run the next line
  # totalReads <- as.integer(as.vector(tumour_stat[,'DP']))
  mutCount =  as.integer(unlist(lapply(strsplit(as.vector(tumour_stat[,'AD']),','),'[[',2)))
  # wtCount = totalReads - mutCount
  wtCount = as.integer(unlist(lapply(strsplit(as.vector(tumour_stat[,'AD']),','),'[[',1)))
  counts_table = mutwt2allelecounts(counts.alt=mutCount,
                                    counts.ref=wtCount,
                                    allele.alt=as.character(vcfdat$V5),
                                    allele.ref=as.character(vcfdat$V4))
  
  output = data.frame(as.character(vcfdat[,1]), vcfdat[,2], counts_table, rowSums(counts_table))
  colnames(output) = c("#CHR","POS","Count_A","Count_C","Count_G","Count_T","Good_depth")
  
  write.table(output, file=outfile, sep="\t", quote=F, row.names=F)
}

library(dpclust3p)

args = commandArgs(TRUE)
vcfdat = read.table(args[1],sep='\t',comment.char='#', stringsAsFactors=F)
datacol = as.integer(args[2]) + 10
battenberg_subclones_file = toString(args[3])
battenberg_cellularity_file = toString(args[4])
sex = "male"
is.male = ifelse(sex=="male", T, F)
namecol = 9

# Create a temp Battenberg rho_psi file for preprocessing
cellularity = read.table(battenberg_cellularity_file, header=T, stringsAsFactors=F)$cellularity
rho_psi = data.frame(rho=c(NA, cellularity, NA), psi=rep(NA, 3), distance=rep(NA, 3), is.best=rep(NA, 3))
row.names(rho_psi) = c("ASCAT", "FRAC_GENOME", "REF_SEG")
battenberg_rho_psi_file = "temp_rho_psi.txt"
write.table(rho_psi, file=battenberg_rho_psi_file, quote=F, col.names=T, row.names=T, sep="\t")
rm(rho_psi)

# Create loci file
loci_file = "loci.txt"
createLociFile(vcfdat, loci_file, 1,2,4,5)
# Create allelecounts file
allelecounts_file = "alleleCounts.txt"
createAlleleCountsFile(vcfdat, datacol, namecol, allelecounts_file)
dpinput_file = "allDirichletProcessInfo.txt"
# Run preprocessing
runGetDirichletProcessInfo(loci_file=loci_file,
                           allele_frequencies_file=allelecounts_file,
                           cellularity_file=battenberg_rho_psi_file,
                           subclone_file=battenberg_subclones_file,
                           gender=sex,
                           SNP.phase.file="NA",
                           mut.phase.file="NA",
                           output_file=dpinput_file)



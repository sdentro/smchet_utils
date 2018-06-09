args = commandArgs(T)

vcf_file = args[1]
outfile = args[2]
samplename = args[3]

library(VariantAnnotation)

addVcfInfoCol = function(vcf, data, number, type, description, abbreviation) {
	  i = header(vcf)@header$INFO
  metadata(vcf)$header@header$INFO <- rbind(i, S4Vectors::DataFrame(Number=number, Type=type, Description=description, row.names=abbreviation))
    info(vcf)[,abbreviation] <- data #as(data, "CharacterList")
    return(vcf)
}

v = readVcf(vcf_file, "hg19")
#counts = as.vector(geno(v)$DP[,colnames(geno(v)$DP)=="tumor"])
# Allelic depths for the ref and alt alleles in the order listed
counts = do.call(rbind, (geno(v)$AD[,colnames(geno(v)$AD)==samplename]))
v = addVcfInfoCol(v, counts[,2], 1, "Integer", "Alt allele supporting reads", "t_alt_count")
v = addVcfInfoCol(v, counts[,1], 1, "Integer", "Ref allele supporting reads", "t_ref_count")
writeVcf(v, filename=outfile)

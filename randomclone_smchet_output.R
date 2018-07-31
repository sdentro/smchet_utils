args = commandArgs(T)
smchet_utils_libpath = args[1]
battenberg_subclones_file = args[2] # T4.battenberg.txt
vcf_file = args[3] # T4.mutect.vcf
assignments_file = args[4] # T4_mutation_assignments.txt
subcl_struct_file = args[5] # T4_subclonal_structure.txt

source(file.path(smchet_utils_libpath, "utils.R"))

vcf = read.table(vcf_file, header=F, comment.char="#", stringsAsFactors=F)
assign = read.table(assignments_file, header=T, stringsAsFactors=F)
chrpos=paste0(vcf$V1, "_", vcf$V2)
chrpos_assign=paste0(assign$chr, "_", assign$pos)
missing = setdiff(chrpos, chrpos_assign)
smchet_asign = rep(0, nrow(vcf))
smchet_asign[chrpos %in% chrpos_assign] = assign$cluster
smchet_asign[is.na(smchet_asign)] = 0

# Add cluster for not assigned mutations
struct = read.table(subcl_struct_file, header=T, stringsAsFactors=F)
num_clusters = nrow(struct)
struct = rbind(struct, data.frame(cluster=0, n_ssms=sum(smchet_asign==0 | is.na(smchet_asign)), proportion=0, ccf=0))
struct = struct[with(struct, order(proportion, decreasing=T)),]

# rename to start counting at 1
struct$clusterid = 1:nrow(struct)
smchet_asign_reindex = smchet_asign
for (i in 1:nrow(struct)) {
	smchet_asign_reindex[smchet_asign==struct$cluster[i]] = struct$clusterid[i]
}

# Build the co-clustering matrix
# co.clustering = get.snv.coassignment.matrix(mut.assignment.type, dataset, iter, burn.in)
no.muts = length(smchet_asign_reindex)
# cellularity = max(optima)
co.clustering = array(0,c(no.muts,no.muts))
for(c in struct$clusterid){
  indices = which(smchet_asign_reindex==c)
  co.clustering[indices,indices] = 1
}
diag(co.clustering) = 1

writeChallengeOutput(battenberg_subclones_file, num_clusters, struct[, c(5,2,3)], smchet_asign_reindex, co.clustering)

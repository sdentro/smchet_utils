args = commandArgs(T)
subclones_file = args[1]

calculate_bb_total_cn = function(bb) {
	  return((bb$nMaj1_A+bb$nMin1_A)*bb$frac1_A + ifelse(!is.na(bb$frac2_A), (bb$nMaj2_A+bb$nMin2_A)*bb$frac2_A, 0))
}

calc_ploidy = function(bb) {
	bb$len = bb$endpos/1000-bb$startpos/1000
	bb$total_cn = calculate_bb_total_cn(bb)
	ploidy = sum(bb$total_cn*bb$len) / sum(bb$len)
	return(ploidy)
}

dat = read.table(subclones_file, header=T, stringsAsFactors=F)
cat(paste0(calc_ploidy(dat), "\n"))

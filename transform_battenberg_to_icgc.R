args = commandArgs(T)
infile_subclones = args[1]
outfile = args[2]

# Example: 
# Rscript make_bb_rounded.R 0c7aca3f-e006-4de3-afc2-20b4f727d4fd_subclones.txt 0c7aca3f-e006-4de3-afc2-20b4f727d4fd_rho_and_psi.txt 0c7aca3f-e006-4de3-afc2-20b4f727d4fd_battenberg_rounded_subclonal_cn.txt

dat = read.table(infile_subclones, header=T, stringsAsFactors=F)
dat_new = dat[,c("chr", "startpos", "endpos")]
colnames(dat_new)[1] = "chromosome"
colnames(dat_new)[2] = "start"
colnames(dat_new)[3] = "end"
dat_new$total_cn = round((dat$nMaj1_A+dat$nMin1_A)*dat$frac1_A + ifelse(!is.na(dat$frac2_A), (dat$nMaj2_A+dat$nMin2_A)*dat$frac2_A, 0))
dat_new$major_cn = round(dat$nMaj1_A*dat$frac1_A) + ifelse(!is.na(dat$frac2_A), round(dat$nMaj2_A*dat$frac2_A), 0)
dat_new$minor_cn = round(dat$nMin1_A*dat$frac1_A) + ifelse(!is.na(dat$frac2_A), round(dat$nMin2_A*dat$frac2_A), 0)
dat_new$total_cn = dat_new$major_cn + dat_new$minor_cn
dat_new$ccf = 1

write.table(dat_new, file=outfile, row.names=F, sep="\t", quote=F)

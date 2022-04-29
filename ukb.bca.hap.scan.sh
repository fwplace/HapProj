#*****************************************************
# UKBB Retrospective Analysis
# Genome-wide scan under a fixed window size
#*****************************************************
ssdir=/scratch_space/fwang2/hap

## Hap-GWAS with w=50
mkdir $ssdir/retrospective/hap.w50.s1.f001.caseOnly
for chrid in {1..22}; do
bsub -P SJLIFEWGS -q standard -J w50s1.chr${chrid} \
     -n 20 -R "rusage[mem=4000] span[hosts=1]" \
     -oo $ssdir/retrospective/hap.w50.s1.f001.caseOnly/chr${chrid}.out \
     -eo $ssdir/retrospective/hap.w50.s1.f001.caseOnly/chr${chrid}.err \
     Rscript ~/hap/scripts/gdsCox_retro_wrap2.R \
             ~/hap/data/ukb.bca.hap.gds/ukb.bca.hap.chr${chrid}.gds \
             50 1 hap 0.001 0.501 20 $ssdir/retrospective/hap.w50.s1.f001.caseOnly 500
done

## Hap-GWAS with w=30
mkdir $ssdir/retrospective/hap.w30.s1.f001.caseOnly
for chrid in {1..22}; do
bsub -P SJLIFEWGS -q large_mem -J w30s1.chr${chrid} \
     -n 20 -R "rusage[mem=4000] span[hosts=1]" \
     -oo $ssdir/retrospective/hap.w30.s1.f001.caseOnly/chr${chrid}.out \
     -eo $ssdir/retrospective/hap.w30.s1.f001.caseOnly/chr${chrid}.err \
     Rscript ~/hap/scripts/gdsCox_retro_wrap2.R \
             ~/hap/data/ukb.bca.hap.gds/ukb.bca.hap.chr${chrid}.gds \
            30 1 hap 0.001 0.501 20 $ssdir/retrospective/hap.w30.s1.f001.caseOnly 500
done

## Hap-GWAS with w=20
mkdir $ssdir/retrospective/hap.w20.s1.f001.caseOnly
for chrid in {1..22}; do
bsub -P SJLIFEWGS -q standard -J w20s1.chr${chrid} \
     -n 20 -R "rusage[mem=4000] span[hosts=1]" \
     -oo $ssdir/retrospective/hap.w20.s1.f001.caseOnly/chr${chrid}.out \
     -eo $ssdir/retrospective/hap.w20.s1.f001.caseOnly/chr${chrid}.err \
     Rscript ~/hap/scripts/gdsCox_retro_wrap2.R \
             ~/hap/data/ukb.bca.hap.gds/ukb.bca.hap.chr${chrid}.gds \
            20 1 hap 0.001 0.501 20 $ssdir/retrospective/hap.w20.s1.f001.caseOnly 500
done

## Hap-GWAS with w=10
mkdir $ssdir/retrospective/hap.w10.s1.f001.caseOnly
for chrid in {1..22}; do
bsub -P SJLIFEWGS -q standard -J w10s1.chr${chrid} \
     -n 20 -R "rusage[mem=4000] span[hosts=1]" \
     -oo $ssdir/retrospective/hap.w10.s1.f001.caseOnly/chr${chrid}.out \
     -eo $ssdir/retrospective/hap.w10.s1.f001.caseOnly/chr${chrid}.err \
     Rscript ~/hap/scripts/gdsCox_retro_wrap2.R \
             ~/hap/data/ukb.bca.hap.gds/ukb.bca.hap.chr${chrid}.gds \
            10 1 hap 0.001 0.501 20 $ssdir/retrospective/hap.w10.s1.f001.caseOnly 500
done

## Hap-GWAS with w=5
mkdir $ssdir/retrospective/hap.w5.s1.f001.caseOnly
for chrid in {1..22}; do
bsub -P SJLIFEWGS -q standard -J w5s1.chr${chrid} \
     -n 20 -R "rusage[mem=4000] span[hosts=1]" \
     -oo $ssdir/retrospective/hap.w5.s1.f001.caseOnly/chr${chrid}.out \
     -eo $ssdir/retrospective/hap.w5.s1.f001.caseOnly/chr${chrid}.err \
     Rscript ~/hap/scripts/gdsCox_retro_wrap2.R \
             ~/hap/data/ukb.bca.hap.gds/ukb.bca.hap.chr${chrid}.gds \
            5 1 hap 0.001 0.501 20 $ssdir/retrospective/hap.w5.s1.f001.caseOnly 500
done

## Hap-GWAS with w=100
mkdir $ssdir/retrospective/hap.w100.s1.f001.caseOnly
for chrid in {1..22}; do
bsub -P SJLIFEWGS -q standard -J w100s1.chr${chrid} \
     -n 20 -R "rusage[mem=4000] span[hosts=1]" \
     -oo $ssdir/retrospective/hap.w100.s1.f001.caseOnly/chr${chrid}.out \
     -eo $ssdir/retrospective/hap.w100.s1.f001.caseOnly/chr${chrid}.err \
     Rscript ~/hap/scripts/gdsCox_retro_wrap2.R \
             ~/hap/data/ukb.bca.hap.gds/ukb.bca.hap.chr${chrid}.gds \
            100 1 hap 0.001 0.501 20 $ssdir/retrospective/hap.w100.s1.f001.caseOnly 500
done

## Hap-GWAS with w=250
mkdir $ssdir/retrospective/hap.w250.s1.f001.caseOnly
for chrid in {1..22}; do
bsub -P SJLIFEWGS -q standard -J w250s1.chr${chrid} \
     -n 20 -R "rusage[mem=4000] span[hosts=1]" \
     -oo $ssdir/retrospective/hap.w250.s1.f001.caseOnly/chr${chrid}.out \
     -eo $ssdir/retrospective/hap.w250.s1.f001.caseOnly/chr${chrid}.err \
     Rscript ~/hap/scripts/gdsCox_retro_wrap2.R \
             ~/hap/data/ukb.bca.hap.gds/ukb.bca.hap.chr${chrid}.gds \
             250 1 hap 0.001 0.501 20 $ssdir/retrospective/hap.w250.s1.f001.caseOnly 500
done

## Hap-GWAS with w=500
mkdir $ssdir/retrospective/hap.w500.s1.f001.caseOnly
for chrid in {1..22}; do
bsub -P SJLIFEWGS -q standard -J w500s1.chr${chrid} \
     -n 20 -R "rusage[mem=4000] span[hosts=1]" \
     -oo $ssdir/retrospective/hap.w500.s1.f001.caseOnly/chr${chrid}.out \
     -eo $ssdir/retrospective/hap.w500.s1.f001.caseOnly/chr${chrid}.err \
     Rscript ~/hap/scripts/gdsCox_retro_wrap2.R \
             ~/hap/data/ukb.bca.hap.gds/ukb.bca.hap.chr${chrid}.gds \
             500 1 hap 0.001 0.501 20 $ssdir/retrospective/hap.w500.s1.f001.caseOnly 500
done

## Summary: initial list iof haplotypes with 10^-7 
nano ~/ProjHap/breast/scripts/hap_retro_sum.R
tmpPath=/scratch_space/fwang2/breast/ukb
for w in 5 10 20 30 50 100 250 500; do
bsub -P ProjHap -q priority -J w${w} -n 1 -R "rusage[mem=8000] span[hosts=1]" \
     -oo $tmpPath/retro/sum.w${w}.out -eo $tmpPath/retro/sum.w${w}.err \
     Rscript ~/ProjHap/breast/scripts/hap_retro_sum.R \
     $tmpPath/retro/hap.w${w}.s1.f001.caseOnly rdata 1e-7 $tmpPath/retro/w${w}.best $tmpPath/retro/w${w}.p7
done

#*****************************************************
# UKBB Prospective Analysis
# Select haplotypes with P<5e-8 for replication
#*****************************************************
nano ~/ProjHap/breast/scripts/hap_prosp_rep.R
tmpPath=/scratch_space/fwang2/breast/ukb
mkdir $tmpPath/prosp
cd $tmpPath/prosp

for w in 5 10 20 30 50 100 250 500; do
mkdir $tmpPath/prosp/w${w}
bsub -P ProjHap -q standard -J w${w} -n 10 -R "rusage[mem=4000] span[hosts=1]" \
     -oo $tmpPath/prosp/w${w}/rep.out -eo $tmpPath/prosp/w${w}/rep.err \
     Rscript ~/ProjHap/breast/scripts/hap_prosp_rep.R 10 \
             $tmpPath/retro/w${w}.p7.rdata $tmpPath/prosp/w${w}/prosp.p5e8.res.data 
done

## combine all windows
R 
df <- NULL
for(w in c(5, 10, 20, 30, 50, 100, 250, 500)){
   load(paste0('w',w,'/prosp.p5e8.res.data'))
   df <- rbind(df, data.frame(run=w, res))
   rm(res)
   print(w)
}
df <- df[, -c(11,14)]
colnames(df)[9:13] <- paste0('retro.',colnames(df)[9:13])
rownames(df) <- 1:nrow(df)
df$rep <- as.numeric(df$pros.p<0.01&sign(log2(df$retro.hr))==sign(log2(df$pros.hr)))
bca.pros.comb <- df
save(bca.pros.comb, file='bca.pros.comb.rdata')
assoc <- bca.pros.comb[bca.pros.comb$rep==1, c('run', 'hid', 'retro.p')]
colnames(assoc) <- c('run','SNP','P')
write.table(assoc, file='bca.pros.rep.assoc', sep='\t',col.names=T, row.names=F, quote=F)

## LD clumping
nano ~/ProjHap/breast/scripts/hap_LDclump.R
tmpPath=/scratch_space/fwang2/breast/ukb
### haplotype data in PLINK file
bsub -P ProjHap -q priority -J LD -n 8 -R "rusage[mem=8000] span[hosts=1]" \
     -oo $tmpPath/prosp/geno.out -eo $tmpPath/prosp/geno.err \
     Rscript ~/ProjHap/breast/scripts/hap_LDclump.R 8 100 $tmpPath/prosp/bca.pros.rep.rdata $tmpPath/prosp

cat << EOF > merge.list
set2.ped set2.map
set3.ped set3.map
set4.ped set4.map
set5.ped set5.map
EOF
plink --file set1 --merge-list merge.list --make-bed --out bca.pros.rep
plink --bfile bca.pros.rep  --clump bca.pros.rep.assoc --clump-p1 5e-8 --clump-p2 5e-8 --clump-r2 0.10 --clump-kb 500 --out bca.pros.rep

for w in 5 10 20 30 50 100 250;do
awk -v ww=${w} 'NR==1 {print; next} NR>1 {if($1==ww) print $1,$2,$3}' bca.pros.rep.assoc > w${w}/w${w}.rep.assoc
plink --bfile bca.pros.rep  --clump w${w}/w${w}.rep.assoc --clump-p1 5e-8 --clump-p2 5e-8 --clump-r2 0.10 --clump-kb 500 --out w${w}/w${w}.rep
done

## get LD and sorted by prospective values
cd /scratch_space/fwang2/breast/ukb
load('bca.pros.comb.rdata')
top <- bca.pros.comb[bca.pros.comb$rep==1, ]
source('~/ProjHap/scripts/ld_find.R')
df <- list(); idx <- 1
for(w in c(5, 10, 20, 30, 50, 100, 250)){
    df[[idx]] <- ld_find(top[top$run==w, ], paste0('w',w,'/w',w,'.rep.clumped'))
    idx <- idx +1
}
df <- do.call('rbind', df)
df2 <- ld_find(df, 'bca.pros.rep.clumped')
df$LD <- df2$ld
df$Lead <- df2$lead
## order LD by prospective p values
df <- df[order(df$pros.p, df$chr, df$start), ]
tmp1 <- unique(df$LD)
tmp2 <- tmp3 <- NULL
for(i in 1:length(tmp1)){
    tmp <- df[df$LD==tmp1[i], ]
    tmp2[i] <- nrow(tmp)
    tmp3[i] <- paste(tmp$chr[1],":",min(tmp$start), '-',max(tmp$end),sep='')    
    print(i)
}
tt <- data.frame(LD=tmp1, bid=1:length(tmp1), LDsize=tmp2, LDregion=tmp3)
tt <- tt[match(df$LD, tt$LD), ]
df$bid <- as.numeric(tt$bid)
df$LDsize <- as.numeric(tt$LDsize)
df$LDregion <- as.vector(tt$LDregion)
df <- df[order(df$bid, df$start), ]
bca.pros.rep <- df[, c(23:25,22,1,19,20,2:17)]
save(bca.pros.rep, file='bca.pros.rep.rdata')

cd /Volumes/clusterhome/fwang2/ProjHap/breast/ukb
w=500
fn <- paste0('retro/w',w, '.best.rdata')
load(fn)
assign('res', get(gsub('.rdata','',basename(fn))))
dim(res)
sum(res$tot)
mean((res$end-res$start)/1e3)

#*****************************************************
# SNP-GWAS
#*****************************************************
nano ~/ProjHap/breast/scripts/ukb_snp.test.R
tmpPath=/scratch_space/fwang2/breast/ukb
mkdir $tmpPath/typed
for i in {1..22};do
bsub -P ProjHap -q standard -J ch${i} -n 10 -R "rusage[mem=4000] span[hosts=1]" \
     -oo $tmpPath/typed/chr${i}.out -eo $tmpPath/typed/chr${i}.err \
     Rscript ~/ProjHap/breast/scripts/ukb_snp.test.R ${i} 10 $tmpPath/typed/snp_fit_chr${i}.rdata
done

## summary
load('/research/projects/yasuigrp/EpiGenetics/common/yasuigrp_data/biobank/gds/ukb.bca.hap/ukb.bca.hap.qc.anno.rdata')
anno <- ukb.bca.hap.qc.anno[, c(1:6,8:14)]
colnames(anno)[11:13] <- c('aachange','cytoband','gwascatalog')
anno$id <- gsub(':','.', anno$id)
anno$id <- gsub('_','.', anno$id)
anno$gene <- gsub('NONE;','', anno$gene)
anno$gene <- gsub(':NONE','', anno$gene)
tt <- list()
for(i in c(1:22)){
   load(paste0('typed/snp_fit_chr',i,'.rdata'))
   t1 <- do.call('rbind', fit)
   t2 <- anno[anno$chr==i, ] 
   t1 <- t1[match(t2$vid, t1$vid), ]
   tt[[i]] <- data.frame(t2, t1[, -1])
   rm(fit)
   print(i)
}
tt <- do.call('rbind', tt)
ukb.bca.typed <- tt
save(ukb.bca.typed, file='ukb.bca.typed.rdata')

#*****************************************************
# UKBB combined analysis of replicated haplotypes
#*****************************************************
cd /scratch_space/fwang2/breast/ukb
nano ~/ProjHap/breast/scripts/ukb_hap.test.R
tmpPath=/scratch_space/fwang2/breast/ukb
bsub -P ProjHap -q priority -J comb -n 10 -R "rusage[mem=4000] span[hosts=1]" \
     -oo ukb.comb.out -eo ukb.comb.err \
     Rscript ~/ProjHap/breast/scripts/ukb_hap.test.R \
         $tmpPath/prosp/bca.pros.rep.rdata 10 $tmpPath/prosp/bca.pros.rep.comb.rdata

## re-calculate LD using combined p
R
load('bca.pros.rep.comb.rdata')
res <- do.call('rbind', fit)
assoc <- res[, c('run', 'hid', 'comb.cox')]
colnames(assoc) <- c('run','SNP','P')
write.table(assoc, file='bca.pros.comb.assoc', sep='\t',col.names=T, row.names=F, quote=F)
## PLINK LD-clumping
plink --bfile bca.pros.rep  --clump bca.pros.comb.assoc --clump-p1 5e-8 --clump-p2 5e-8 --clump-r2 0.10 --clump-kb 500 --out bca.pros.comb
## summary
load('bca.pros.rep.comb.rdata')
res <- do.call('rbind', fit)
ld_find <- function(dat, fn.ld){
    tt <- read.table(file=fn.ld, header=T)[, c(1:5,12)]
    tt <- tt[order(tt$CHR, tt$BP), ]
    sp2 <- list()
    for(j in 1:nrow(tt)){
          if(tt$SP2[j]=='NONE') sp2[[j]] <- tt$SNP[j]
          if(tt$SP2[j]!='NONE'){
               tmp <- strsplit(tt$SP2[[j]], '[,]')[[1]]
               tmp <- as.vector(unlist(sapply(tmp, function(x){strsplit(x,split='[()]')[[1]][1]})))
               sp2[[j]] <- c(tt$SNP[j], tmp)
          }
     }
     ld <- data.frame(ld=rep(1:nrow(tt), lapply(sp2,length)),
                      sp1=rep(tt$SNP,lapply(sp2,length)),sp2=as.vector(unlist(sp2)))
     ld <- ld[match(dat$hid, ld$sp2), ]
     dat$LD <- ld$ld
     tt <- dat[order(dat$comb.cox), ]
     tt <- tt[!duplicated(tt$LD), ]
     dat$Lead <- as.numeric(dat$hid%in%tt$hid)
    return(dat)
}
res2 <- ld_find(res, 'bca.pros.comb.clumped')
## combined p does not change the LD of 20 replicated loci
bca.pros.rep.comb <- res2[, 1:32]
save(bca.pros.rep.comb, file='bca.pros.rep.comb.rdata')

## copy data to home directory
cp -R ukb ~/ProjHap/breast

#*****************************************************
# dbGap validation
#*****************************************************
cd /Volumes/clusterhome/fwang2/ProjHap/breast/ukb
R
## SNP members of each haplotype
load('prosp/bca.pros.rep.comb.rdata')
load('ukb.bca.typed.rdata')
load('/Volumes/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/ukb/qced/dbgap.ukb.qced.snp.rdata')
hapdat  <- bca.pros.rep.comb
snpdat  <- ukb.bca.typed
hap2snp <- list()
for(i in 1:nrow(hapdat)){
  tt <- as.numeric(strsplit(hapdat$block[i], '[_]')[[1]])
  dat <- snpdat[snpdat$chr==tt[1], ]
  dat <- dat[dat$vid%in%tt[2]:tt[3], ]
  dat <- data.frame(hid=hapdat$hid[i], dat)
  dat$idx <- 1:nrow(dat)
  rownames(dat) <- dat$idx
  dat$aachange <- as.vector(unlist(sapply(dat$aachange,function(x){strsplit(x,'[,]')[[1]][1]})))
  dat$allele <- as.integer(strsplit(hapdat$hap[i],'')[[1]][-1])
  ## imputation QC
  vv <- dbgap.ukb.qced.snp[dbgap.ukb.qced.snp$chr==dat$chr[1], ]
  vv <- vv[match(dat$id, vv$id.hg19), ]
  dat$imp.rsq  <- vv$rsq
  dat$imp.miss <- vv$miss
  dat$imp.hwe  <- vv$hwe
  dat$imp.maf  <- vv$maf
  hap2snp[[i]] <- dat
  print(i)
  rm(tt, dat, vv)
}
ukb.bca.rep20 <- list(hapdat=hapdat, hap2snp=hap2snp)
save(ukb.bca.rep20, file='ukb.bca.rep20.rdata')

## SNP reduction
tmpPath=/scratch_space/fwang2/breast/ukb
mkdir $tmpPath/reduction
for i in {1..436};do
bsub -P ProjHap -q standard -J h${i} -n 1 -R "rusage[mem=4000] span[hosts=1]" \
     -oo $tmpPath/reduction/h${i}.out -eo $tmpPath/reduction/h${i}.err \
     Rscript ~/ProjHap/breast/scripts/ukb_hap_reduction.R ${i} 0.80 \
             ~/ProjHap/breast/ukb/ukb.bca.rep20.rdata $tmpPath/reduction/hap${i}.rdata 
done

./R
df <- list()
for(i in 1:436){
   load(paste0('/scratch_space/fwang2/breast/ukb/reduction/hap',i,'.rdata'))
   df[[i]] <- reduction.obj
   rm(reduction.obj)
   print(i)
}
ukb.bca.rep20.reduction.rsq80 <- df
save(ukb.bca.rep20.reduction.rsq80, file='ukb.bca.rep20.reduction.rsq80.rdata')

## validation 1: SNP reduction with HDS-Max
cd /home/fwang2/ProjHap/breast/ukb
mkdir validation

load('ukb.bca.rep20.reduction.rsq80.rdata')
sets <- lapply(ukb.bca.rep20.reduction.rsq80, function(x){x$snpset})
save(sets, file='validation/topmed.rsq80.snpsets.rdata')

tmpPath=/scratch_space/fwang2/breast/ukb
mkdir $tmpPath/validation
mkdir $tmpPath/validation/max0

run='max0'
mkdir $tmpPath/validation/${run}
for i in {1..436};do
bsub -P PrjHap -q standard -J h${i} -n 4 -R "rusage[mem=4000] span[hosts=1]" \
     -oo $tmpPath/validation/${run}/h${i}.out -eo $tmpPath/validation/${run}/h${i}.err \
      Rscript ~/ProjHap/breast/scripts/haptest_reduction_exempt0.R ~/ProjHap/breast/ukb/validation/pheno.overall.txt \
      /research/projects/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/ukb/qced/dbgap.topmedr2.ukb.xxx.qced.hg19.gds \
      ~/ProjHap/breast/ukb/validation/topmed.rsq80.snpsets.rdata max ${i} 4 $tmpPath/validation/${run}/h${i}.rdata
done

run='sum0'
mkdir $tmpPath/validation/${run}
for i in {1..436};do
bsub -P PrjHap -q standard -J h${i} -n 4 -R "rusage[mem=4000] span[hosts=1]" \
     -oo $tmpPath/validation/${run}/h${i}.out -eo $tmpPath/validation/${run}/h${i}.err \
      Rscript ~/ProjHap/breast/scripts/haptest_reduction_exempt0.R ~/ProjHap/breast/ukb/validation/pheno.overall.txt \
      /research/projects/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/ukb/qced/dbgap.topmedr2.ukb.xxx.qced.hg19.gds \
      ~/ProjHap/breast/ukb/validation/topmed.rsq80.snpsets.rdata sum ${i} 4 $tmpPath/validation/${run}/h${i}.rdata
done

## validation 2: SNP reduction with HDS-Sum2
load('ukb.bca.rep20.rdata')
load('ukb.bca.rep20.reduction.rsq80.rdata')
hdat <- ukb.bca.rep20$hapdat
hidx <- which(hdat$retro.afcase<0.01)
## 193 rare haplotypes
ukb.bca.rare13 <- list(hdat=hdat[hidx, ], hap2snp=ukb.bca.rep20$hap2snp[hidx])
sets <- lapply(ukb.bca.rep20.reduction.rsq80[hidx], function(x){x$snpset})
save(sets, file='validation/topmed.rsq80.rare13.snpsets.rdata')
save(ukb.bca.rare13, file='ukb.bca.rare13.rdata')

run='sum1'
mkdir $tmpPath/validation/${run}
for i in {1..193};do
bsub -P PrjHap -q standard -J h${i} -n 4 -R "rusage[mem=4000] span[hosts=1]" \
     -oo $tmpPath/validation/${run}/h${i}.out -eo $tmpPath/validation/${run}/h${i}.err \
      Rscript ~/ProjHap/breast/scripts/haptest_reduction_exempt1.R ~/ProjHap/breast/ukb/validation/pheno.overall.txt \
      /research/projects/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/ukb/qced/dbgap.topmedr2.ukb.xxx.qced.hg19.gds \
      ~/ProjHap/breast/ukb/validation/topmed.rsq80.rare13.snpsets.rdata sum ${i} 4 $tmpPath/validation/${run}/rare_h${i}.rdata
done

run='sum2'
mkdir $tmpPath/validation/${run}
for i in {1..193};do
bsub -P PrjHap -q standard -J h${i} -n 4 -R "rusage[mem=4000] span[hosts=1]" \
     -oo $tmpPath/validation/${run}/h${i}.out -eo $tmpPath/validation/${run}/h${i}.err \
      Rscript ~/ProjHap/breast/scripts/haptest_reduction_exempt2.R ~/ProjHap/breast/ukb/validation/pheno.overall.txt \
      /research/projects/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/ukb/qced/dbgap.topmedr2.ukb.xxx.qced.hg19.gds \
      ~/ProjHap/breast/ukb/validation/topmed.rsq80.rare13.snpsets.rdata sum ${i} 4 $tmpPath/validation/${run}/rare_h${i}.rdata
done

## Summary
tmpPath=/scratch_space/fwang2/breast/ukb
cd $tmpPath/validation

nano  ~/ProjHap/breast/scripts/combine_res.R
args <- commandArgs(trailingOnly = TRUE)
get_sum <- function(respath, reduction){
  fn <- list.files(respath, 'rdata', full.names=T)
  tmp <- gsub('rare_h','', gsub('.rdata','',basename(fn)))
  tmp <- gsub('h','', tmp)
  tmp <- as.numeric(tmp)
  fn <- fn[order(as.numeric(tmp))]
  df <- list()
  for(i in 1:length(fn)){
     load(fn[i])
     if(reduction){
           if(class(res[[1]])=='list') tmp1 <- lapply(res, function(x){do.call('rbind', x)})
           if(class(res[[1]])=='data.frame') tmp1 <- res
           idx <- which(sapply(tmp1, function(x){class(x)=='data.frame'}[1]))
           tmp2 <- as.numeric(sapply(tmp1[idx], nrow))
           tt <- data.frame(msplit=rep(c(NA, 5,10,20,40)[idx], tmp2), do.call('rbind', tmp1[idx]))
           tt$lgp <- ifelse(tt$OR>1, 1, -1)*-log10(apply(tt[, c(12:13)],1,min))
           df[[i]] <- tt
     }else{
           if(class(res)=='list') tt <- do.call('rbind', res)
           if(class(res)=='data.frame') tt <- res
           tt$lgp <- ifelse(tt$OR>1, 1, -1)*-log10(apply(tt[, c(11:12)],1,min))
           df[[i]] <- tt
     }
     rm(res, tt)
     print(i)
  }
return(df) 
}
assign(as.character(args[3]), get_sum(as.character(args[1]), args[2]))
save(list=as.character(args[3]), file=paste0(as.character(args[3]),'.rdata'))

Rscript ~/ProjHap/breast/scripts/combine_res.R max0 TRUE  topmed.rsq80.max0
Rscript ~/ProjHap/breast/scripts/combine_res.R sum0 TRUE  topmed.rsq80.sum0
Rscript ~/ProjHap/breast/scripts/combine_res.R sum1 TRUE  topmed.rsq80.sum1
Rscript ~/ProjHap/breast/scripts/combine_res.R sum2 TRUE  topmed.rsq80.sum2

cp * ~/ProjHap/breast/ukb/validation

cd /Volumes/clusterhome/fwang2/ProjHap/breast/ukb
load('ukb.bca.rep20.rdata')
hdat <- ukb.bca.rep20$hapdat
load('ukb.bca.rep20.reduction.rsq80.rdata')
load('validation/topmed.rsq80.max0.rdata')
load('validation/topmed.rsq80.sum0.rdata')
topmed.rsq80.val1 <- list()
for(i in 1:436){
    t1 <- topmed.rsq80.max0[[i]]
    t2 <- topmed.rsq80.sum0[[i]]
    t1 <- data.frame(snpexempt=0, t1)
    t2 <- data.frame(snpexempt=0, t2)
    topmed.rsq80.val1[[i]] <- list(sum=ukb.bca.rep20.reduction.rsq80[[i]]$sum, res=rbind(t1,t2))
  print(i)
}
save(topmed.rsq80.val1, file='validation/topmed.rsq80.val1.rdata')

load('ukb.bca.rep20.rdata')
hdat <- ukb.bca.rep20$hapdat
load('ukb.bca.rep20.reduction.rsq80.rdata')
load('validation/topmed.rsq80.sum0.rdata')
load('validation/topmed.rsq80.sum1.rdata')
load('validation/topmed.rsq80.sum2.rdata')
hidx <- which(hdat$retro.afcase<0.01)
ukb.bca.rep20.reduction.rsq80 <- ukb.bca.rep20.reduction.rsq80[hidx]
topmed.rsq80.sum0 <- topmed.rsq80.sum0[hidx]
topmed.rsq80.val2 <- list()
for(i in 1:193){
    t1 <- topmed.rsq80.sum0[[i]]
    t2 <- topmed.rsq80.sum1[[i]]
    t3 <- topmed.rsq80.sum2[[i]]
    t1 <- data.frame(snpexempt=0, t1)
    t2 <- data.frame(snpexempt=1, t2)
    t3 <- data.frame(snpexempt=2, t3)
    topmed.rsq80.val2[[i]] <- list(sum=ukb.bca.rep20.reduction.rsq80[[i]]$sum, res=rbind(t1, t2, t3))
  print(i)
}
save(topmed.rsq80.val2, file='validation/topmed.rsq80.val2.rdata')

load('ukb.bca.rep20.rdata')
load('ukb.bca.rep20.reduction.rsq80.rdata')
hdat <- ukb.bca.rep20$hapdat
hidx <- which(hdat$retro.afcase<0.01)
ukb.bca.rare13.reduction <- ukb.bca.rep20.reduction.rsq80[hidx]
save(ukb.bca.rare13.reduction, file='ukb.bca.rare13.reduction.rdata')


#*****************************************************
# Tables and Plots
#*****************************************************
cd ~/fwork/ProjHap/breast/ukb
mkdir res
R
load('ukb.bca.rep20.rdata')
hdat <- ukb.bca.rep20$hapdat
df <- do.call('rbind', lapply(ukb.bca.rep20$hap2snp, function(x){x[which.min(x$comb.cox), ]}))
df <- df[, c(1,7,2:6,8,10:14,28,27,29:32)]
colnames(df) <- gsub('comb.', 'snp.', colnames(df))
com <- data.frame(hdat[, -c(6:8,17,23,29)], df[, -c(1,4,14,15)])
write.table(colnames(com), file='res/top436_header.txt',sep='\n', col.names=F, row.names=F, quote=F)
cols <- as.vector(unlist(read.table(file='res/top436_header.txt',header=F)))
com <- com[, cols]
for(i in c(5,6,10,11,15,16,25)){
    com[, i] <- 1e4*com[, i]
}
ukb.bca.top436 <- com
save(ukb.bca.top436, file='res/ukb.bca.top436.rdata')
write.csv(com, file='res/ukb.bca.top436.csv')

## summary of Validation
load('ukb.bca.top436.rdata')
load('validation/topmed.rsq80.val1.rdata')
ukb.bca.top436$imp <- sapply(topmed.rsq80.val1, function(x){x$sum$reduction[1]})
hidx <- which(ukb.bca.top436$Lead==1)
df <- ukb.bca.top436[hidx, c(1:3, 33:36,38)]

tt <- NULL
for(i in hidx){
     t1 <- topmed.rsq80.val1[[i]]$sum
     t2 <- topmed.rsq80.val1[[i]]$res
     t1 <- t1[!is.na(t1$msplit)&t1$msplit==40, ]
     t2 <- t2[!is.na(t2$msplit)&t2$msplit==40&t2$type=='max'&t2$hcut==0.500, ]
     tt <- rbind(tt, data.frame(t1[, c(1,3,7,9,8,11,12)], t2[, -c(2)]))
}
r1 <- data.frame(val='max0 wo c', df, tt)

tt <- NULL
for(i in hidx){
     t1 <- topmed.rsq80.val1[[i]]$sum
     t2 <- topmed.rsq80.val1[[i]]$res
     t1 <- t1[!is.na(t1$msplit)&t1$msplit==40, ]
     t2 <- t2[!is.na(t2$msplit)&t2$msplit==40&t2$type=='sum', ]
     t2 <- t2[which.max(t2$lgp), ]
     tt <- rbind(tt, data.frame(t1[, c(1,3,7,9,8,11,12)], t2[, -c(2)]))
}
r2 <- data.frame(val='sum0 wi c', df, tt)

tt <- rbind(r1,r2)
tt$ukb.afconc <- 1e4*tt$ukb.afconc
tt$ukb.afcase <- 1e4*tt$ukb.afcase
write.csv(tt, file='res/dbgap_val1.csv')

load('validation/topmed.rsq80.val2.rdata')
bids <- c(5,6,9,11:20)
df <- ukb.bca.top436[ukb.bca.top436$bid%in%bids, c(1:3,32:36,38)]
tt <- NULL
for(i in 1:193){
     t1 <- topmed.rsq80.val2[[i]]$sum
     t2 <- topmed.rsq80.val2[[i]]$res
     t2 <- t2[!is.na(t2$msplit)&t2$msplit%in%c(10,40)&t2$type=='sum'&!is.na(t2$snpidx)&t2$snpexempt==1, ]
     t2 <- t2[which.max(t2$lgp), ]
     t1 <- t1[!is.na(t1$msplit)&t1$msplit==t2$msplit, ]
     tt <- rbind(tt, data.frame(t1[, c(1,3,7,9,8,11,12)], t2[, -c(2)]))
}
r1 <- data.frame(val='sum1 wi c', df, tt)

tt <- NULL
for(i in 1:193){
     t1 <- topmed.rsq80.val2[[i]]$sum
     t2 <- topmed.rsq80.val2[[i]]$res
     t2 <- t2[!is.na(t2$msplit)&t2$msplit%in%c(10,40)&t2$type=='sum'&!is.na(t2$snpidx)&t2$snpexempt==2, ]
     t2 <- t2[which.max(t2$lgp), ]
     t1 <- t1[!is.na(t1$msplit)&t1$msplit==t2$msplit, ]
     tt <- rbind(tt, data.frame(t1[, c(1,3,7,9,8,11,12)], t2[, -c(2)]))
}
r2 <- data.frame(val='sum2 wi c', df, tt)
tt <- rbind(r1,r2)
tt$ukb.afconc <- 1e4*tt$ukb.afconc
tt$ukb.afcase <- 1e4*tt$ukb.afcase
write.csv(tt, file='res/dbgap_val2.csv')
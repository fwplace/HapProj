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

## Summary: initial list iof haplotypes with 5x10^-8 
nano ~/ProjHap/breast/scripts/hap_retro_sum.R
tmpPath=/scratch_space/fwang2/breast/ukb
for w in 5 10 20 30 50 100 250 500; do
bsub -P ProjHap -q priority -J w${w} -n 1 -R "rusage[mem=8000] span[hosts=1]" \
     -oo $tmpPath/retro/sum.w${w}.out -eo $tmpPath/retro/sum.w${w}.err \
     Rscript ~/ProjHap/breast/scripts/hap_retro_sum.R \
     $tmpPath/retro/hap.w${w}.s1.f001.caseOnly rdata 5e-8 $tmpPath/retro/w${w}.best $tmpPath/retro/w${w}.p5e8
done

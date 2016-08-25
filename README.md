# admix_from_vcf
Run admixture on illumina VCF files using python and 1000 genomes phase 3 PED and MAP files. 

## Requires:
- Plink 1.9
- VCFtools 0.1.13
- SnpEff 4.2 (build 2015-12-05)
- Gzipped individual vcf files, illumina format
 
## 1000 Genomes Phase 3 Reference Files 
Downolad the following files and update their location on your computer, as needed: 
Files were created as outlined here:
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/admixture_files/README.admixture_20141217
PED file:
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/admixture_files/ALL.wgs.phase3_shapeit2_filtered.20141217.maf0.05.ped
MAP file:
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/admixture_files/ALL.wgs.phase3_shapeit2_filtered.20141217.maf0.05.map
Panel file:
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
References:
LD Pruning: http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#prune
Pre-Formatted Files are located at:
/users/selasady/shared/admix
run_admix()
- creates reference marker file set
- converts VCF of interest to plink bed/bim/fam format
- runs admixture
- saves results to an admixture_results.txt file
"""

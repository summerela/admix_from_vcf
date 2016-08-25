#!/usr/bin/env python

"""
Requires:
- Plink 1.9
- VCFtools 0.1.13
- SnpEff 4.2 (build 2015-12-05)
- Gzipped individual vcf files, illumina format

#### 1000 Genomes Phase 3 Reference Files #####
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

####################
### Script Setup ###
####################

# import modules
import os
import subprocess as sp
import pandas as pd


# disable extraneous pandas warning
pd.options.mode.chained_assignment = None

class run_admix(object):

    # program paths
    tool_path = '/users/selasady/my_titan_itmi/tools/'
    plink_path = "{}plink1.9/plink".format(tool_path)
    java_path = '/tools/java/jdk1.7/bin/java -jar -Xmx16g'
    gatk_path = "{}GenomeAnalysisTK.jar".format(tool_path)
    admixture_path = "{}admixture_linux-1.3.0/admixture".format(tool_path)
    vcf_verify = "{}snpEff/scripts/vcfBareBones.pl".format(tool_path)
    bcftools = "{}bcftools/bcftools".format(tool_path)
    vcftools = "{}vcftools_0.1.13/bin/vcftools".format(tool_path)

    def __init__(self, vcf_dir, ref_dir, ref_base, panel_file, population='super', num_cores=15):
        """Returns an admix object with related parameters and files."""
        self.vcf_dir = vcf_dir
        self.ref_dir = ref_dir
        self.ref_base = ref_base
        self.panel_file = panel_file
        self.population = population
        self.num_cores = num_cores

    #################################################################################
    ## setup basic methods for working with files and calling subprocess commands  ##
    #################################################################################
    @staticmethod
    # create function to run bash command with subprocess
    def subprocess_cmd(command, input_dir=os.getcwd()):
        '''
        Run programs in bash via subprocess
        :param command: command string as would be run on the command line
        :param input_dir: optional directory to run command in, default cwd
        :return: runs bash command
        '''
        print ("Running \n {}".format(command))
        ps = sp.Popen(command, shell=True, stdout=sp.PIPE, stderr=sp.PIPE, cwd=input_dir)
        try:
            print ps.communicate()
        except sp.CalledProcessError as e:
            print e

    @staticmethod
    def get_file_base(in_file, num_to_remove):
        basename = str('.'.join(in_file.split('.')[:-int(num_to_remove)]) if '.' in in_file else in_file)
        return basename

    @staticmethod
    def create_output_dir(self, input_dir):
        if not os.path.exists(input_dir):
            os.makedirs(input_dir)

    ################################################################
    ### Create reference marker bed/bim/fam and supporting files ###
    ################################################################

    # create panel text file with subject id's, grouped by population
    def create_sorter_and_pop(self):
        '''
        Creates sorted popuplation files for reference data set for super and subpopulations
        :param panel_file: patht to tab-delimited .panel file in format [subject_id | super_pop | sub_pop  | gender]
        :return sorted_ref.txt a two column tsv ( sample | pop ) to be fed to plink
        for sorting by populations and super_pop.txt and sub_pop.txt
        '''
        out_panel = "{}/sorted_panel.txt".format(self.ref_dir)

        try:
            # read in 1000g admix panel file
            panel = pd.read_csv(self.panel_file, sep='\t', index_col=False, header=0, usecols=[0,1,2,3])
            # sort the panel by super and then sub pop for ease of downstream analysis
            panel.sort_values(by=['super_pop', 'pop'], axis=0, inplace=True)
            # create super_pop.txt file for super populations
            super_out = panel[['sample', 'super_pop']]
            print ("Creating super_pop.txt...")
            super_out.to_csv("{}/super_pop.txt".format(self.ref_dir), index=False,
                             names=['sample', 'super_pop'], sep='\t', header=False)
            # create sub_pop.txt for subpopulations
            sub_out = panel[['sample', 'pop']]
            print ("Creating sub_pop.txt...")
            sub_out.to_csv("{}/sub_pop.txt".format(self.ref_dir), index=False, names=['sample', 'sub_pop'],
                           sep='\t', header=False)
            # create sorted_ref.txt for sorting marker tsv by pop in format family_id | subject_id
            sorter_out = panel[['sample', 'sample']]
            print("Saving sorted panel file as {}".format(out_panel))
            sorter_out.to_csv("{}/sorted_panel.txt".format(self.ref_dir), index=False, sep='\t', header=False)
        except Exception as e:
            raise SystemExit("Exception: {}".format(e))

    def filter_ped_bed(self):
        '''
        Filter 1000g ped file (requires map file in same dir) by linkage disequilibrium and MAF
        :return: Marker lists written to filtered_maf_ld.prune.in and filtered_maf_ld.prune.out
        to be bed into the bed/bim/bam creation for keeping only this filtered marker set
        '''
        filter_ped_cmd = "{plink} --file {ref_dir}/{base} --out filtered_maf_ld --maf 0.05 --max-maf .49 \
            --indep-pairwise 50 5 0.5".format(ref_dir=self.ref_dir, plink=self.plink_path, base=self.ref_base)
        self.subprocess_cmd(filter_ped_cmd)

    def ped2bed(self):
        '''
        Convert the ped and map files to bed/bim/fam required for admixture, using the sorted panel to sort the output
        and the pruned file filtered by MAF and LD using filter_ped_bed() to keep only filtered markers
        :return: bed, bim, fam files from ped/map files
        '''
        pruned_file = '{ref_dir}/filtered_maf_ld.prune.in'.format(ref_dir=self.ref_dir)
        ped2bed_cmd = '{plink} --file {ref}/{base} --extract {pruned} --make-bed --indiv-sort file {ref}/sorted_panel.txt \
                --out {ref}/1000g_marker'.format(plink=run_admix.plink_path, ref=self.ref_dir, base=self.ref_base,
                                                   pruned=pruned_file)
        self.subprocess_cmd(ped2bed_cmd)

    def make_snp_list(self):
        '''
        Takes MAF and LD filtered bim file and creates markfer file called snplist.txt
        :return: snplist.txt of filtered marker variants
        '''
        print ("Creating marker file at {}/snplist.txt".format(self.ref_dir))
        bim_file = "{ref_dir}/1000g_marker.bim".format(ref_dir=self.ref_dir)
        # pull out chrom and position from bim file
        ref_bim = pd.read_csv(bim_file, sep='\t', header=None, usecols=[0, 3], names=['chrom', 'pos'])
        # sort by chrom and pos
        ref_bim.sort_values(by=['chrom', 'pos'], axis=0, inplace=True)
        # add chr prefix to match with itmi vcf format
        ref_bim['chrom'] = 'chr' + ref_bim['chrom'].astype(str)
        # save to a file
        ref_bim.to_csv("{}/snplist.txt".format(self.ref_dir), header=False, index=False, sep='\t')

    ################################################
    ### Pull out matching markers from VCF files ###
    ################################################
    def subset_vcf(self, input_vcf):
        '''
        Pulls out marker positions from gzipped vcf file and strip extraneous info
        :param input_vcf: gzipped vcf file in vcf_dir
        :return: gzipped vcf file containing only chrom, pos, ref and alt at marker positions
        '''
        out_base = self.get_file_base(input_vcf, 2)
        snp_file = "{}/snplist.txt".format(self.ref_dir)
        subset_cmd = r'''{vcftools} --gzvcf {vcf_dir}/{input_vcf} --positions {snps} --recode --stdout | {barebones} |
        gzip > {out_dir}/{out_vcf}_trimmed.vcf.gz'''\
            .format(vcftools=self.vcftools, vcf_dir=self.vcf_dir, input_vcf=input_vcf,
                    snps=snp_file, barebones=self.vcf_verify, out_dir=self.vcf_dir,
                    out_vcf=out_base)
        self.subprocess_cmd(subset_cmd)

    def vcf_to_bed(self, input_vcf):
        '''
        convert vcf to bed/bim/fam with plink 1.9 with following options:
        - double-id = Set both FIDs and IIDs to the VCF/BCF sample ID
        - biallelic-only strict = keep only biallelic variants
        - allow-no-sex = allow subjects with no gender info
        - set-missing-var-ids = keeps variants with no rsid by assigning them one
        :param input_vcf: vcf processed with subset_vcf function above
        :return: subject.genome_trimmed.bed/bim/fam output to vcf_dir
        '''
        bed_out = self.get_file_base(input_vcf, 2)
        vcf2plink_cmd = "nohup {plink} --vcf {vcf_dir}/{file} --double-id --biallelic-only strict \
        --allow-no-sex --geno 0.999 --set-missing-var-ids @:#[b37]\$1,\$2 --make-bed --out {vcf_dir}/{bed}"\
            .format(plink=self.plink_path, vcf_dir=self.vcf_dir, file=input_vcf, bed=bed_out)
        print ("Converting {} from vcf to bed/bim/fam").format(input_vcf)
        self.subprocess_cmd(vcf2plink_cmd)

    ######################E###################
    ### Merge subject beds with Marker bed ###
    ##########################################

    def add_subject_to_sorter(self, input_bed):
        '''
        Add subject id to sorter file to populations in order while creating
        merged bed files
        :param input_bed: input bed file to merge with marker
        :return: sampleID_sorted.txt sorted file with reference and individual subject ids
        '''
        base_name = self.get_file_base(input_bed, 3)
        old_sorter = pd.read_csv("{}/sorted_panel.txt".format(self.ref_dir), sep='\t', header=None,
                                 names=['family', 'sample'])
        new_sorter = old_sorter.append({'family': base_name, 'sample': base_name}, ignore_index=True)
        out_name = "{}/{}_sorted.txt".format(self.vcf_dir, base_name)
        print("Saving new sorter panel to {}".format(out_name))
        new_sorter.to_csv(out_name, sep='\t', index=None)

    def find_beds(self, input_vcf):
        '''
        Locate bed files in a dir and ensure they have matching bim/fam
        :param input_dir: dir where bed files to merge are located
        :return: list of bed files to merge
        '''
        bedList = []
        base_name = self.get_file_base(input_vcf, 2)
        bedList.append(str(self.vcf_dir + '/' + base_name + '.bed'))
        bedList.append(str(self.vcf_dir + '/' + base_name + '.bim'))
        bedList.append(str(self.vcf_dir + '/' + base_name + '.fam'))
        bedList.sort()
        return bedList

    def merge_beds(self, input_vcf):
        '''
        Generate command to merge bed files using plink then run with subprocess
        :param input_vcf:
        :return: merged bed/bim/fam file of 1000g marker plus input bed files
        '''
        print ("\n Merging bed file(s) with marker file... \n")
        base_name = self.get_file_base(input_vcf, 1)
        bed_base = self.get_file_base(input_vcf, 3)
        beds = self.find_beds(input_vcf)
        bed_list = ' '.join(beds)
        merge_cmd = "{plink} --bfile {ref_dir}/1000g_marker --bmerge {bed_list}  --indiv-sort  file \
        {out_dir}/{bed_base}_sorted.txt --memory 40000 \
        --out {out_dir}/{bed_base}_merged".format \
            (plink=self.plink_path, ref_dir=self.ref_dir, bed_list=bed_list, ref=self.ref_dir,
             bed_base=bed_base, vcf=base_name, out_dir=self.vcf_dir)
        self.subprocess_cmd(merge_cmd)



    #######################
    ###  run admixture  ###
    #######################

    def make_pop_file(self, input_bed):
        '''
        Create .pop file with 5 known major populations from 1000g and '-' for unknown for each merged bed
        Assumes ref_dir contains tsv files in format subject_id | pop named super_pop.txt and sub_pop.txt
        these are created with create_sorter_and_pop() function
        :param input_vcf: input_bed file to grab base file name
        :param kval: specify 5 for super population and 26 for subpopulations
        :return .pop file for merged bed/bim/fam set as input for supervised admixture
        '''
        # read in _merged.fam file
        base_name = self.get_file_base(input_bed, 1)
        in_fam = "{}/{}.fam".format(self.vcf_dir, base_name)
        merged_df = pd.read_csv(in_fam, sep='\t', usecols=[1], names=['sample'])
        if self.population == 'super':
            # merge .fam file with pop file to get population for each known subject
            super_in = "{}/super_pop.txt".format(self.ref_dir)
            super_map = pd.read_csv(super_in, sep='\t', header=None, names=['sample', 'super_pop'])
            super_out = pd.merge(merged_df, super_map, how='left', on='sample')
            super_out.drop('sample', inplace=True, axis=1)
            super_file = "{}/{}.pop".format(self.vcf_dir, base_name)
            print ("Saving pop file as {}".format(super_file))
            super_out.to_csv(super_file, sep='\t', na_rep='-', header=False, index=False)
        elif self.population == 'sub':
            # create sub_pop file
            sub_in = "{}/sub_pop.txt".format(self.ref_dir)
            sub_map = pd.read_csv(sub_in, sep='\t', header=None, names=['sample', 'sub_pop'])
            sub_out = pd.merge(merged_df, sub_map, how='left', on='sample')
            sub_out.drop('sample', inplace=True, axis=1)
            sub_file = "{}/{}.pop".format(self.vcf_dir, base_name)
            print ("Saving pop file as {}".format(sub_file))
            sub_out.to_csv(sub_file, sep='\t', na_rep='-', header=False, index=False)
        else:
            print ("Please set population to 'super' or 'sub' when instantiating admix class.")

    def run_admix_program(self, merged_bed):
        '''
        run admixture on input merged bed files
        :input: bed file merged from reference and individual bed/bim/fam files
        requires supporting pop file created with make_pop_file() function
        :return: p and q files for determining admixture
        '''
        if self.population == 'super':
            print("Running admixture with super populations for {}".format(merged_bed))
            super_cmd = "nohup {} -j4 {} --supervised 5".format(self.admixture_path, merged_bed)
            self.subprocess_cmd(super_cmd)
        elif self.population == 'sub':
            print("Running admixture with sub-populations for {}".format(merged_bed))
            sub_cmd = "nohup {} -j4 {} --supervised 26".format(self.admixture_path, merged_bed)
            self.subprocess_cmd(sub_cmd)

    ######################
    ## Process Results ###
    ######################

    def get_results(self, q_file):
        # read admixture results in
        q_base = self.get_file_base(q_file, 2)
        q_file_in = pd.read_csv("{}/{}".format(admix.vcf_dir, q_file), header=None, sep=' ')
        pop_in = "{}/{}.pop".format(admix.vcf_dir, q_base)
        # read in pop file
        pop_file = pd.read_csv(pop_in, header=None, sep=' ', names=['pop'])
        # read in fam file
        fam_file_in = "{}/{}.fam".format(admix.vcf_dir, admix.get_file_base(q_file, 2))
        fam_file = pd.read_csv(fam_file_in, header=None, sep='\t', usecols=[1], names=['sample'])
        frames = [fam_file, pop_file, q_file_in]
        out_frame = pd.concat(frames, ignore_index=True, keys=None, axis=1)
        pop_list = list(out_frame[1][(out_frame[1] != '-')].unique())
        head_list = ['sample', 'known_pop']
        col_list = head_list + pop_list
        out_frame.columns = col_list
        results_df = out_frame[(out_frame['known_pop'] == '-')]
        results_df.drop(['known_pop'], inplace=True, axis=1, errors='ignore')
        results_df.to_csv("{}/admix_results.txt".format(admix.vcf_dir, q_base), header=True, sep='\t', index=False,
                      mode='a')


#################################################
if __name__ == '__main__':

    # update these options with each run
    admix = run_admix(vcf_dir='/users/selasady/shared/admix/vcf_files/', \
                      panel_file='/users/selasady/shared/admix/integrated_call_samples_v3.20130502.ALL.panel', \
                      ref_dir='/users/selasady/shared/admix', \
                      ref_base='ALL.wgs.phase3_shapeit2_filtered.20141217.maf0.05',
                      population='super')

#######################################

    # TODO comment this out if marker file already created
    # create files for 1000g markers
    admix.create_sorter_and_pop()
    admix.filter_ped_bed()
    admix.ped2bed()
    admix.make_snp_list()

    # create genome_trimmed.vcf
    for file in os.listdir(admix.vcf_dir):
        if file.endswith('.vcf.gz'):
            admix.subset_vcf(file)

    # create genome_trimmed.bed/bim/fam
    # add subject id to sorter file
    # create sorted bed/bim/fam with ref + sample_id
    for file in os.listdir(admix.vcf_dir):
        if file.endswith('_trimmed.vcf.gz'):
            admix.vcf_to_bed(file)
            admix.add_subject_to_sorter(file)
            admix.merge_beds(file)

    # running merged file through admix
    for file in os.listdir(admix.vcf_dir):
        if file.endswith('_merged.bed'):
            admix.make_pop_file(file)
            admix.run_admix_program(file)

    # process results of admix
    for file in os.listdir(admix.vcf_dir):
        if file.endswith('.Q'):
            admix.get_results(file)




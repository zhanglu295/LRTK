# Description of LRTK-SEQ

`LRTK-SEQ` (Linked Reads Toolkit for resequencing) is an all-in-one package and designed to analyze linked reads sequencing data from 10X Chromium system. We implemented several functions to staisfify the requirements of the resequencing analysis for human genome.

## Requirements
### Software

`LRTK-SEQ` was implemented by [`python3`](https://www.python.org/downloads/source/) and three external packages are needed to be pre-installed: [`numpy`](http://www.numpy.org/), [`pysam`](http://pysam.readthedocs.io/en/latest/index.html), [`matplotlib`](https://matplotlib.org/). LRTK-SEQ also applied several available programes [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [`Picard`](http://broadinstitute.github.io/picard/), [`GATK`](https://software.broadinstitute.org/gatk/download/auth?package=GATK), [`BWA`](https://github.com/lh3/bwa), [`SAMtools`](http://www.htslib.org/), [`sambamba`](http://lomereiter.github.io/sambamba/), [`sbt`](http://www.scala-sbt.org/), [`fgbio`](https://github.com/fulcrumgenomics/fgbio), [`HapCut2`](https://github.com/vibansal/HapCUT2) and [`NAIBR`](https://github.com/raphael-group/NAIBR). If some software are already installed in your system, please make a softlink of the executable program to the software directory that specified by `LRTK-SEQ.py Config -s` 

### Dataset
reference genome(fasta) <br>
barcode whitelist(fasta) <br>
https://software.broadinstitute.org/gatk/download/bundle

`LRTK-SEQ` also provides an auxiliary programe `preinstall.sh` to preinstall all the requisite automatically. Please make sure your system is connected to the internet. All the software and dataset are downloaded and installed to `./bin` and `./dataset`, respectively.

## Donwload `LRTK-SEQ` to local machine

Please use `--recursive` to clone submodules simultaneously.

git clone --recursive https://github.com/zhanglu295/LRTK

## Quick start
Step1: install LRTK-SEQ

`git clone --recursive https://github.com/zhanglu295/LRTK`

Step2: install softare and dataset

`sh preinstall.sh`

Step3: generating configure files

`python3 LRTK-SEQ.py Config -o outputdir`

Step4: Linked read preprocessing and alignment

`python3 LRTK-SEQ.py Basicall -i fqlist.txt -c outputdir/config/Basic.config -o outputdir`

Step5: variant calling and phasing

`python3 LRTK-SEQ.py Reseqall -i basic.bam.list -c outputdir/config/Reseq.config -o outputdir`

## Basic usage

python3 LRTK-SEQ.py \<command\> [options]

### Command:
    Config		Generate configuration files
    Basicall	Run the whole pipeline of data preprocessing and alignment, including PRE, MFQ, ALN, MARK, BQSR, STAT and MERGE 
    Reseqall	Run the whole pipeline of variant calling and phasing, including Varcall(SNV/InDel, based on individual), FBV(SNV/InDel, based on family) and Phasing
    Clean		Delete temporary files

    Basic		Run customized steps in `Basicall`
    Reseq		Run customized steps in `Reseqall`
    
#### Config
    python3 LRTK-SEQ.py Config [options]
    
    Basic options:
    -o --outputdir, output directory.
    
    Advanced options:
    -s --softwarepath, string, software directory [default: {abs_path(LRTK-SEQ.py)}/bin]
    -d --datasetpath, string, dataset directory [default: {abs_path(LRTK-SEQ.py)}/dataset]
    -b --bed, string, bed of genome regions without Ns [default: {abs_path(LRTK-SEQ.py)}/dataset/GATK_bundle/nonN.bed]
    
`Config` would generate two configuration files: 1. Basic.config for the components in `Basicall` 2. Reseq.config for the components in `Reseqall`, and intervals file used in `Reseqall` <br>

#### Basicall
    python3 LRTK-SEQ.py Basicall [0-8] [options]
    
    eg. python3 LRTK-SEQ.py Basicall [options]  (run all steps in Basicall step by step)
    eg. python3 LRTK-SEQ.py Basicall 1-3 [options] (run steps 1, 2, and 3 in Basicall)
    eg. python3 LRTK-SEQ.py Basicall 4- [options] (run steps 4, 5, 6, and 7 in Basicall)
    step 1: PRE     generate clean fastq files and correct barcode error
    step 2: MFQ     merge fq based on library or sample
    step 3: ALN     alignment
    step 4: MARK    barcode-aware PCR duplicates removal
    step 5: BQSR    base quality score recalibration
    step 6: STAT    calculate QC statistics based on sample library, including Cf, Cr, MuFL, NFP etc.
    step 7: MERGE   merge all bam files from the same sample

    Warnings: The input file must change according to the starting steps!
    Warnings: The input file must change according to the starting steps!!
    Warnings: The input file must change according to the starting steps!!!
    For instance, if "Basicall" starts with step 3, the input file must change to be the input file of "ALN"

    Considering the enormous fastq files generated by whole genome deep sequecning, the original fastq files would be splitted into smaller ones to imporve efficiency, in the steps `PRE`, `MFQ` and `ALN`.

    Basic options:
    -i --input, string, input file for fastq information (4 or 5 columns:1.Sample ID; 2.Library ID; 3. Serial number of the library; 4.Path to the fastqs [5. Path to the fastq2]), or others(according to the starting step).
    -o --outputdir, string, the path to output
    -c --config, string, configuration file [default: outputdir/config/Basic.config]

    Advanced options:
    -M --mergefq, int [0, 1 or 2. default: 0]. 0: do nothing, MFQ would not run; 1: merge fq files belong to the same library; 2: merge fq files belong to the same sample. If '-M' > 0, all following analysis would be based on the merged file. Other numbers might cause unpredictable error resuls.
    -m --minlen, int, the minimum length (bp) of molecule to be considered [default: 500]
    -p --parallel, int, the number CPUs can be used in parallel [default: 1]
    -z --splitsize, int, the amount of reads of splited fastq, reads num = split_size / 8 [default: 60000000 lines, 7500000 read-pairs, compressed file size: ~300M]

`Basicall` includes seven submodules: `PRE`, `MFQ`, `ALN`, `MARK`, `BQSR`, `STAT` and `MERGE`, which would be carried out step by step. Users can also choose any steps to run independently by using `‘Basic’`. <br><br>
`Basicall` would generate the processed fastq files, recalibrated bam files, and filtered (based on molecule info) bam files which would be used for variation calling.

eg. <br>
python3 LRTK-SEQ.py Basicall -i fqlist.txt -o outputdir -c outputdir/config/Basic.config <br>

###### Basic:
    python3 LRTK-SEQ.py Basic <command> [options]
    
    PRE       generate clean fastq files and correct barcode error
    MFQ       merge fq based on library or sample
    ALN       alignment
    MARK      barcode-aware PCR duplicates removal
    BQSR      base quality score recalibration
    STAT      calculate QC statistics based on sample library, including Cf, Cr, MuFL, NFP etc
    MERGE     merge all bam files from the same sample

###### PRE:
    python3 LRTK-SEQ.py Basic PRE [options]

    Basic options:
    -i --input, string, input file for fastq information (4 or 5 columns:1.Sample ID; 2.Library ID; 3. Serial number of the library; 4.Path to the fastqs [5. Path to the fastq2])
    -o --outputdir, string, the path to output
    -c --config, string, configuration file [default: outputdir/config/Basic.config]

    Advanced options:
    -p --parallel, int, the number CPUs can be used in parallel [default: 4]
    -z --splitsize, int, the amount of reads of splited fastq, reads num = split_size / 8 [default: 60000000 lines, 7500000 read-pairs, compressed file size: ~300M]

`PRE` is similar as what `basic` in [`Long ranger`](https://support.10xgenomics.com/genome-exome/software/pipelines/latest/what-is-long-ranger) do, including barcode error correction, barcode white-listing, attaching barcodes to reads. The description of output files are written in outputdir/Result_list/Basic_PRE_result.txt. Meanwhile, it also generates the input file for the next step `MFQ` or `ALN`: outputdir/Result_list/Basic_MFQ_input.txt or outputdir/Result_list/Basic_ALN_input.txt <br><br>

eg. <br>
python3 LRTK-SEQ.py Basic PRE -i fqlist.txt -o outputdir -c outputdir/config/Basic.config -p 2 -z 80000000<br>

###### MFQ:
    python3 LRTK-SEQ.py Basic MFQ [options]

    Basic options:
    -i --input, string, input file for fastq information (4 columns:1.Sample ID; 2.Library ID; 3. fastq list; 4. prefix of new barcode, including 4 letters made up by "A", "T", "C", or "G"), generated by 'PRE'
    -o --outputdir, string, the path to output
    -c --config, string, configuration file [default: outputdir/config/Basic.config]

    Advanced option:
    -M --mergefq, int [0, 1 or 2. default: 0]. 0: do nothing, MFQ would not run; 1: merge fq files belong to the same library; 2: merge fq files belong to the same sample. If '-M' > 0, all following analysis would be based on the merged file. Other numbers might cause unpredictable error resuls.
    -p --parallele, int, the number CPUs can be used in parallel [default: 4] 

`Basic MFQ` would merge fastq files belong to the sample library(-M 1) or sample(-M 2), or do nothing(-M 0). The barcode of each reads would be re-named in the processing of merging. The description of output files are written in outputdir/Result_list/Basic_MFQ_result.txt. Meanwhile, it also generates the input file for the next step `ALN`: outputdir/Result_list/Basic_ALN_input.txt <br><br>

eg. <br>
python3 LRTK-SEQ.py Basic PRE -i outputdir/Result_list/Basic_MFQ_input.txt -o outputdir -c outputdir/config/Basic.config -p 2 -M 2<br>

###### ALN:
    python3 LRTK-SEQ.py Basic ALN [options]
    
    Basic options:
    -i --input, string, input file for fastq information (4 columns:1.Sample ID; 2.Library ID; 3. Serial number of the library; 4.fastqs file list generated in `PRE` or `MFQ`), generated by 'PRE' or 'MFQ'
    -o --outputdir, string, the path to output
    -c --config, string, configuration file [default: outputdir/config/Basic.config]

    Advanced options
    -p --parallel, the number CPUs can be used in parallel [default: 4] 

`ALN` is similar as what `align` in [`Long ranger`](https://support.10xgenomics.com/genome-exome/software/pipelines/latest/what-is-long-ranger) do, read alignment by `bwa`. The description of output files are written in outputdir/Result_list/Basic_ALN_result.txt. Meanwhile, it also generates the input file for the next step `MARK`: outputdir/Result_list/Basic_MARK_input.txt. <br><br>

eg. <br>
python3 LRTK-SEQ.py Basic ALN -i outputdir/Result_list/Basic_ALN_input.txt -o outputdir -c outputdir/config/Basic.config <br>

###### MARK
    python3 LRTK-SEQ.py Basic MARK [options]
    
    Basic options:
    -i --input, string, input file for BAM information (3 columns: 1.Sample Id; 2.Library Id; 3.Path to BAM, generated by `ALN`)
    -o --outputdir, string, the path to output
    -c --config, string, configuration file [default: outputdir/config/Basic.config]

`MARK` would merge all the bams from the same library and perform barcode-aware PCR duplicates removal using [`picard`]. The description of output files are written in outputdir/Result_list/Basic_MARK_result.txt. Meanwhile, it also generates the input file for the next step `BQSR`: outputdir/Result_list/Basic_BQSR_input.txt. <br><br>

e.g. <br>
python3 LRTK-SEQ.py Basic MARK -i outputdir/Result_list/Basic_MARK_input.txt -o outputdir -c outputdir/config/Basic.config <br>

###### BQSR
    python3 LRTK-SEQ.py Basic BQSR [options]
    
    Basic options:
    -i --input, string, input file for BAM information (3 columns: 1.Sample Id; 2.Library Id; 3.Path to BAM, generated by `MARK`)
    -o --outputdir, string, the path to output
    -c --config, string, configuration file [default: outputdir/config/Basic.config]

`BQSR` performs base quality score recalibration in [`GATK`](https://software.broadinstitute.org/gatk/).

The description of output files are written in outputdir/Result_list/Basic_BQSR_result.txt. Meanwhile, it also generates the input file for the next step `STAT`: outputdir/Result_list/Basic_STAT_input.txt. <br><br>

e.g.<br>
python3 LRTK-SEQ.py Basic BQSR -i outputdir/Result_list/Basic_BQSR_input.txt -o outputdir -c outputdir/config/Basic.config <br>

###### STAT
    python3 LRTK-SEQ.py Basic STAT [options]
    
    Basic options
    -i --input, string, input file for BAM information (3 columns: 1.Sample Id; 2.Library Id; 3.Path to BAM, generated by `BQSR`)
    -o --outputdir, string, the path to output
    -c --config, string, configuration file [default: outputdir/config/Basic.config]

    Advanced options
    -m --minlen, int, the minimum length (bp) of molecule to be considered [default: 500]

The short reads align to the adjacent locations and share the same barcodes should come from the same long DNA fragment. `LRTK-SEQ` reconstructs these long DNA fragments and writes their information into `fragment.csv`. It contains nine columns: 1. fragment id 2. barcode sequence 3. chromosome 4. start position 5. end position 6. fragment length 7. number of reads 8. total length of all the reads in this fragment 9. depth coverage. 

`LRTK-SEQ` also generates four histograms: 1. Unweighted fragment length distribution 2. Weighted fragment length distribution 3. number of fragments per droplet 4. The distribution of sequencing depth per fragment.

The description of output files are written in outputdir/Result_list/Basic_STAT_result.txt. Meanwhile, it also generates the input file for the next step `MERGE`: outputdir/Result_list/Basic_MERGE_input.txt. <br><br>

e.g. <br>
python3 LRTK-SEQ.py Basic STAT -i outputdir/Result_list/Basic_STAT_input.txt -o outputdir -c outputdir/config/Basic.config <br>

###### MERGE
    python3 LRTK-SEQ.py Basic MERGE [options]
    
    Basic options:
    -i --input, string, input file for BAM information (3 columns: 1.Sample Id; 2.Library Id; 3.Path to BAM, generated by `MARK`)
    -o --outputdir, string, the path to output
    -c --config, string, configuration file [default: outputdir/config/Basic.config]
 
`MERGE` merges all bam files belong to the same sample. The description of output files are written in outputdir/Result_list/Basic_MERGE_result.txt. Meanwhile, it also generates the input file for the next module `Reseq`: outputdir/Result_list/Reseq_Varcall_input.txt. <br><br>

e.g. <br>
python3 LRTK-SEQ.py Basic MERGE -i outputdir/Result_list/Basic_MERGE_input.txt -o outputdir -c outputdir/config/Basic.config <br>

#### Reseqall
    python3 LRTK-SEQ.py Reseqall [options]
    
    Basic options:
    -i --input, string, input file for BAM information (2 columns: 1. Sample id 2. Path to BAM, generated by `Basicall`)
    -o --outputdir, string, the path to output
    -c --config, string, configuration file [default: outputdir/config/Reseq.config]
    -f --familyinfo, string, family information, including 2 columns: 1. family id, 2. sample id

`Reseqall` includes three submodules: `Varcall`, `Famcall` and `Phasing`, which would be carried out step by step automatically. Users can also choose any steps to run independently by using `‘Reseq’`.

`Reseqall` would generate unphased VCF files of each individual/family, and phased VCF files of each individual.

eg. <br>
python3 LRTK-SEQ.py Reseqall -i outputdir/Result_list/Reseq_Varcall_input.txt -o outputdir -c outputdir/config/Reseq.config <br>

##### Reseq
    python LRTK-SEQ.py Reseq <command> [options]
    
    Varcall     call SNVs and InDels by GATK, based on individual
    Famcall     call SNVs and InDels by GATK, based on family
    Phasing     phasing variants by HapCUT2

###### Varcall
    python3 LRTK-SEQ.py Reseq Varcall [options]
    
    Basic options:
    -i --input, string, input file for BAM information (2 columns: 1. Sample id 2. Path to BAM), generated by `Basicall`
    -o --outputdir, string, the path to output
    -c --config, string, configuration file [default: outputdir/config/Reseq.config]  

`Varcall` detects SNVs and InDels using GATK4, based on individuals. The output vcf files are written in outputdir/Result_list/Reseq_Varcall_result.txt, and input for `Famcall` and `Phasing` is written in outputdir/Result_list/Reseq_Famcall_input.txt and outputdir/Result_list/Reseq_Phasing_input.txt, respectively. <br><br>

e.g.
python3 LRTK-SEQ.py Reseq Varcall -i outputdir/Result_list/Reseq_Varcall_Phasing_input.txt -o outputdir -c outputdir/config/Reseq.config <br>

###### Famcall
    python3 LRTK-SEQ.py Reseq Famcall [options]

    Basic options:
    -i --input, string, input file for variation calling (2 columns: 1. Sample id 2. gvcf file list), generated by `Varcall`
    -o --outputdir, string, the path to output
    -c --config, string, configuration file [default: outputdir/config/Reseq.config]
    -f --familyinfo, string, family information, including 2 columns: 1. family id, 2. sample id

`Famcall` detects SNVs and InDels using GATK4, based on families. The output vcf files are written in outputdir/Result_list/Reseq_Famcall_result.txt. <br><br>

e.g.
python3 LRTK-SEQ.py Reseq Famcall -i outputdir/Result_list/Reseq_Famcall_input.txt -o outputdir -c outputdir/config/Reseq.config -f family_info.txt <br>

###### Phasing
    python3 LRTK-SEQ.py Reseq Phasing [options]
    
    Basic options:
    -i --input, string, input file for phasing (2 columns: 1. Sample id 2. unphased vcf), generated by `Varcall`
    -o --outputdir, string, the path to output
    -c --config, string, configuration file [default: ./config/Reseq.config]

`Phasing` phases the variants in vcf files by [HapCUT2](https://github.com/vibansal/HapCUT2). The phased vcf is written in outputdir/Result_list/Reseq_Phasing_result.txt .<br><br>

e.g. <br>
python LRTK-SEQ.py Reseq Phasing -i outputdir/Result_list/Reseq_Varcall_input.txt -o outputdir -c outputdir/config/Reseq.config <br>



`Reseqall` includes three submodules: `Varcall`, `SVcall` and `Phasing`, which would be carried out step by step automatically. Users can also choose any steps to run independently by using `‘Reseq’`. output?

eg. <br>
python LRTK-SEQ.py Reseqall -i out_dir/Result_list/Reseq_Varcall_Phasing_input.txt -o out_dir -c ./config/Reseq.config <br>

##### Reseq
    python LRTK-SEQ.py Reseq <command> [options]
    
    Varcall     call SNVs and Indels by GATK
    SVcall      call structure variantions by NAIBR 
    Phasing     phasing variants by HapCUT2

###### Varcall
    python LRTK-SEQ.py Reseq Varcall [options]
    
    Basic options:
    -i --input, the BAM information files from Basicall or ALN
    -o --outputdir, the path to output
    -c --config, configuration file [default: ./config/Reseq.config]
 
    Advanced optionL
    -p --parallel, the number CPUs can be used in parallel [default: 1]
    

`Varcall` detects SNVs and InDels using GATK3. The output vcf files are written in out_dir/Result_list/Reseq_Varcall_result.txt, and input for `Phasing` is written in out_dir/Result_list/Reseq_Phasing_input.txt. <br><br>

e.g.
python LRTK-SEQ.py Reseq Varcall -i out_dir/Result_list/Reseq_Varcall_Phasing_input.txt -o out_dir -c ./config/Reseq.config <br>

###### SVcall
    python LRTK-SEQ.py Reseq SVcall [options]
    
    Basic options:
    -i --input, the BAM information files from Basicall or ALN
    -o --outputdir, the path to output
    -c --config, configuration file [default: ./config/Reseq.config]
 
    Advanced optionL
    -p --parallel, the number CPUs can be used in parallel [default: 1]
    -m --min_mapq: Minimum mapping quality for a read to be included in analysis (default: 40)
    -s --min_sv: Minimum size of a structural variant to be detected (default: lmax, the 95th percentile of the paired-end read insert size distribution)
    -k --min_barcode minimum number of barcode overlaps supporting a candidate NA (default = 3)
    

`SVcall` detects structure variations by NAIBR. The output vcf files are written in out_dir/Result_list/Reseq_SVcall_result.txt. <br><br>

e.g.
python LRTK-SEQ.py Reseq Varcall -i out_dir/Result_list/Reseq_Varcall_Phasing_input.txt -o out_dir -c ./config/Reseq.config <br>

###### Phasing
    python LRTK-SEQ.py Reseq Phasing [options]
    
    Basic options:
    -i --input, the BAM information files from Basicall or ALN
    -v --vcf, unphased vcf file generated by Varcall
    -o --outputdir, the path to output
    -c --config, configuration file [default: ./config/Reseq.config]

`Phasing` phases the variants in vcf files by [HapCUT2](https://github.com/vibansal/HapCUT2). The phased vcf is written in out_dir/XX.vcf.<br><br>

e.g. <br>
python LRTK-SEQ.py Reseq Phasing -i out_dir/Result_list/Reseq_Varcall_Phasing_input.txt -v out_dir/Result_list/Reseq_Phasing_input.txt -o out_dir -c ./config/Reseq.config <br>

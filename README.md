# Description of LRTK-SEQ

`LRTK-SEQ` (Linked Reads Toolkit for resequencing) is an all-in-one package and designed to analyze linked reads sequencing from 10X Chromium System. We implemented several essential submodules to staisfify the requirements of resequencing data analysis for human genome.

## Requirements

LRTK-SEQ was implemented by [python3](https://www.python.org/downloads/source/) and three external packages are needed to be pre-installed: [numpy](http://www.numpy.org/), [pysam](http://pysam.readthedocs.io/en/latest/index.html), [matplotlib](https://matplotlib.org/)

## Donwload `LRTK-SEQ` to local machine

Please use `--recursive` to clone essential software simultaneously.

git clone --recursive https://github.com/zhanglu295/LRTK


## Basic usage

python LRTK-SEQ.py \<command\> [options]

### Command:
    Config		Generate configuration file
    Basicall	Run the whole pipeline of basic data preparation, including CFQ_ALN, MARK, BQSR, STAT, MERGE 
    Reseqall	Run the whole pipeline of resequencing, including Varcall, SVcall, Phasing
    Clean		Delete temporary files

    Basic		Run customized steps for basic data preparation
    Reseq		Run customized steps for resequencing

## Quick start
install LRTK-SEQ
`git clone --recursive https://github.com/zhanglu295/LRTK`
install softare and data
`sh preinstall.sh`


    
#### Config
    python LRTK-SEQ.py Config [options]
    
    Options:
    -o --outputdir, Output directory. (required)
    -s --softwarepath, Software directory [default: ./bin]
    -d --datasetpath, Dataset directory [default: ./dataset]
    -b --bed, Target regions for whole-exome/Target sequencing, or non-N regions for whole genome sequencing [default: ./dataset/genome/nonN.bed]. If nonN.bed is missing, it would be generated automatically with respect to the reference sequence in ‘./dataset/genome/genome.fa’
    
`Config` would generate two configuration files, one is for `Basicall`/`Basic` (outputdir/config/Basic.config) and the other is for `Reseqall`/`Reseq` (outputdir/config/Reseq.config). Meanwhile, the list of chromosomes (used for parallelization) would be listed in XXX according to dataset/genome/genome.fa.fai. <br>

Required datasets in the `dataset` folder are: 
1. indexed genome fasta files(dataset/genome/genome.fa)
2. indexed barcode fasta files (dataset/barcode/barcode.fa). 

Certainly, any changes of these file or file names are allowded, you just need to guarantee the correspongdings of all your input files in the final configuration file.<br>


#### Basicall
    python LRTK.py Basicall [options]
    
    requisite options:
    -i --input, the input file containing fastq information (The input file contains three columns:1.Sample ID;2.Library ID;3. Path to sample fastqs).
    -o --outputdir, the path to output
    -c --config, configuration file [default: outdir/config/Basic.config]

    alternative options:
    -m --molecule, number. molecule length less than it would be discarded [default: 500]
    -p --parallel, the number of fq pairs that processing parallel. The max amount of invoking CPU would be 4*(-p) [default: 1]
    -N --noBX, generate additional fq file that has BX info or not [default: yes]

`Basicall` including `CFQ_ALN`, `MARK`, `BQSR`, `STAT` and `MERGE`, it would run all these patterns step by step automatically. Moreover, you can also run these patterns step by step yourself, using `‘Basic’` patterns. <br><br>

eg. <br>
python LRTK.py Basicall -i fqlist.txt -o outputdir -c outputdir/config/Basic.config <br>

##### Basic:
    python LRTK.py Basic \<command\> [options]
    
    CFQ_ALN   generate clean fastq files and correct barcode error
    MARK      merge all bam files belong to the same library of each sample, and barcode aware PCR duplication removal (must complete ALN) 
    BQSR      recalibrate base quality scores, using GATK.
    STAT      calculate QC statistics, including Cf, Cr, MuFL, NFP etc. (must complete CFQ and ALN)
    MERGE     merge all bam files belong to the same sample

###### CFQ_ALN
    python LRTK.py Basic CFQ_ALN
    
    -i --input, the input file containing fastq information (The input file contains three columns:1.Sample ID;2.Library ID;3. Path to sample fastqs).
    -o --outputdir, the path to output
    -p --parallel, the number of fq pairs that processing parallel. The max amount of invoking CPU would be 4*(-p) [default: 1]
    -N --noBX, generate additional fq file that has BX info or not [default: yes]
    -c --config, configuration file [default: outdir/config/Basic.config]

`CFQ_ALN` verifies the barcode info and calculates the sequencing quality of original fastq files, and performs alignment using [`bwa`](http://bio-bwa.sourceforge.net/). In order to speed up the scheme, multithreading is implemented, and the parallel number has been fixed as 4. Output files would be listed in outputdir/Result_list/Basic_CFQ_ALN_result.txt. Meanwhile, list of input files for the next step `MARK` would be generate at the end of `CFQ_ALN`: outputdir/Result_list/Basic_MARK_input.txt. <br><br>

eg. <br>
python LRTK.py Basic CFQ_ALN -i fqlist.txt -o outputdir -c outputdir/config/Basic.config <br>
###### MARK
    python LRTK.py Basic MARK

    -i --input, the input file containing the information of bam files (The input file contains two columns:1.Sample Id;2.Library Id;3.Path to bam)
    -o --outputdir, the path to output
    -p --parallele, the number of fq pairs that processing parallel. The max amount of invoking CPU would be 4*(-p) [default: 1]
    -c --config, configuration file [default: outdir/config/Basic.config]

`MARK` merges bam files belonging to the same library(based on the 2nd column of the input file), and mark duplication reads that generated in PCR using [`picard`](http://broadinstitute.github.io/picard/). Output files of `MARK` would be listed in outputdir/Result_list/Basic_MARK_result.txt, and list of input file for the next step `BQSR` would be generated at the end of `MARK`: outputdir/Result_list/Basic_BQSR_input.txt. <br><br>

e.g. <br>
python LRTK.py Basic MARK -i outputdir/Result_list/Basic_MARK_input.txt -o outputdir -c outputdir/config/Basic.config <br>
###### BQSR
    python LRTK.py Basic BQSR

    -i --input, the input file that contains the SAM/BAM files generated by STAT (The input file contains three columns:1.Sample Id;2.Library Id;3.Path to bam)
    -o --outputdir, the path to output
    -p --parallele, the number of fq pairs that processing parallel. The max amount of invoking CPU would be 4*(-p) [default: 1]
    -c --config, configuration file [default: outdir/config/Basic.config]

`BQSR` performs base quality score recalibration using [`GATK`](https://software.broadinstitute.org/gatk/), corresponding datasets can be download in (https://software.broadinstitute.org/gatk/download/bundle). Parallel number has also been fixed as 4. Output files of `BQSR` would be listed in outputdir/Result_list/Basic_BQSR_result.txt, and list of input files for the next step `STAT` would be generated at the end of `BQSR`: outputdir/Result_list/Basic_STAT_input.txt. <br><br>

e.g.<br>
python LRTK.py Basic BQSR -i outputdir/Result_list/Basic_BQSR_input.txt -o outputdir -c outputdir/config/Basic.config <br>

###### STAT
    python LRTK.py Basic STAT

    -i --input, the input file that contains the SAM/BAM files generated by MAK (The input file contains two columns:1.Sample Id;2.Library Id;3.Path to bam)
    -o --outputdir, the path to output
    -p --parallele, the number of fq pairs that processing parallel. The max amount of invoking CPU would be 4*(-p) [default: 1]
    -m --molecule, number. molecule length less than it would be discarded [default: 500]
    -c --config, configuration file [default: outdir/config/Basic.config]

`STAT` calculates molecule information for each library(based on the 2nd column of input file), you can control the smallest molecule with the help of `-m`[default: 500], and molecules with length smaller than it would be discarded. Meanwhile, some important statistics, for instance, CF and CR et al., would be calculated to evaluate the sequencing quality. Moreover, some figures would be generated, such as Cr, weighted and unweighted molecule length distribution, and so on. <br>
Similarly, output files of `STAT` would be listed in outputdir/Result_list/Basic_STAT_result.txt, and list of input files for the next step `MERGE` would be generated at the end of `STAT`: outputdir/Result_list/Basic_MERGE_input.txt. <br><br>

e.g. <br>
python LRTK.py Basic STAT -i outputdir/Result_list/Basic_STAT_input.txt -o outputdir -c outputdir/config/Basic.config <br>
###### MERGE
    python LRTK.py Basic MERGE
    
    -i --input, the input file that contains the SAM/BAM files generated by MAK (The input file contains two columns:1.Sample Id;2.Library Id;3.Path to bam)
    -o --outputdir, the path to output
    -c --config, configuration file [default: outdir/config/Basic.config]

`MERGE` merges bam files belonging to the same sample, output files would be listed in outputdir/Result_list/Basic_MERGE_result.txt, and list of input files for `Reseq` would be generated at the end of `MERGE`: outputdir/Result_list/Reseq_Varcall_Phasing_input.txt. <br><br>

e.g. <br>
python LRTK.py Basic MERGE -i outputdir/Result_list/Basic_MERGE_input.txt -o outputdir -c outputdir/config/Basic.config <br>

#### Reseqall
    python LRTK.py Reseqall [options]
    
    -i --input, the input file that contains the BAM files generated by Basicall or ALN
    -o --outputdir, the path to output
    -L --chrlist, list of chromosome [e.g. chr1, chr2, chrX, default: outdir/config/chrlist.txt]
    -p --parallel, the number of CPU allowed to use [default: 1]
    -c --config, configuration file [default: outdir/config/Reseq.config]

`Reseqall` including `Varcall` and `Phasing`, it would run all these patterns step by step automatically. Moreover, you can also run these patterns step by step yourself, using `‘Reseq’` patterns. <br><br>

eg. <br>
python LRTK.py Reseqall -i outputdir/Result_list/Reseq_Varcall_Phasing_input.txt -o outputdir -c outputdir/config/Reseq.config <br>

##### Reseq
    python LRTK.py Reseq [options]
    
    Varcall     call SNVs and Indels by GATK
    SVcall      call structure variantion by GROC-SVs (developing, do not work)
    Phasing     phasing variants by HapCUT2

###### Varcall
    python LRTK.py Reseq Varcall

    -i --input, the input file that contains the BAM files generated by Basicall or ALN
    -o --outputdir, the output directory path
    -L --chrlist, list of chromosomes [e.g.chr1,chr2,chrX, default: outdir/config/chrlist.txt]
    -p --parallel, the number of CPU allowed to use [default: 1]
    -c --config, configuration file [default: outdir/config/Reseq.config]

`Varcall` calls SNVs and InDels by chromosome based on GATK3, the pipeline may not work if you use GATK4 or other versions. Original vcf files would be listed in outputdir/Result_list/Reseq_Varcall_result.txt, and input files for the next step `Phasing` would be listed in outputdir/Result_list/Reseq_Phasing_input.txt. <br><br>

e.g.
python LRTK.py Reseq Varcall -i outputdir/Result_list/Reseq_Varcall_Phasing_input.txt -o outputdir -L outputdir/config/chrlist.txt -c outputdir/config/Reseq.config <br>

###### Phasing
    python LRTK.py Reseq Phasing

    -i --input, the input path that contains the BAM files generated by Basicall or ALN
    -v --vcf, unphased vcf file generated by Varcall or the other variant callers, both compressed or uncompressed vcf files are allowed
    -o --outputdir, the path to output
    -c --config, configuration file [default: outdir/config/Reseq.config]

`Phasing` phases vcf files using [HapCUT2](https://github.com/vibansal/HapCUT2).<br><br>

e.g. <br>
python LRTK.py Reseq Phasing -i outputdir/Result_list/Reseq_Varcall_Phasing_input.txt -v outputdir/Result_list/Reseq_Phasing_input.txt -o outputdir -c outputdir/config/Reseq.config <br>

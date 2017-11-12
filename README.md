# Description of LRTK-SEQ

`LRTK-SEQ` (Linked Reads Toolkit for resequencing) is an all-in-one package and designed to analyze linked reads sequencing data from 10X Chromium system. We implemented several functions to staisfify the requirements of the resequencing analysis for human genome.

## Requirements
### Software

`LRTK-SEQ` was implemented by [`python3`](https://www.python.org/downloads/source/) and three external packages are needed to be pre-installed: [`numpy`](http://www.numpy.org/), [`pysam`](http://pysam.readthedocs.io/en/latest/index.html), [`matplotlib`](https://matplotlib.org/). LRTK-SEQ also applied several available programes [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [`Picard`](http://broadinstitute.github.io/picard/), [`GATK`](https://software.broadinstitute.org/gatk/download/auth?package=GATK), [`BWA`](https://github.com/lh3/bwa), [`SAMtools`](http://www.htslib.org/), [`sbt`](http://www.scala-sbt.org/), [`fgbio`](https://github.com/fulcrumgenomics/fgbio), [`HapCut2`](https://github.com/vibansal/HapCUT2) and [`NAIBR`](https://github.com/raphael-group/NAIBR). If some software are already installed in your system, please make a softlink of the executable program to the software directory that specified by `LRTK-SEQ.py Config -s` 

### Dataset
reference genome(fasta)
barcode whitelist(fasta)
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

`python LRTK-SEQ.py Config -o ./config`

Step4: Linked read preprocessing and alignment

`python LRTK-SEQ.py Basicall -i fqlist.txt -c ./config/Basic.config -o out_dir`

Step5: variant calling and phasing

`python LRTK-SEQ.py Reseqall -i basic.bam -o ./reseqout -c ./config/Reseq.config`

## Basic usage

python LRTK-SEQ.py \<command\> [options]

### Command:
    Config		Generate configuration files
    Basicall	Run the whole pipeline of data preprocessing and alignment, including CFQ_ALN, MARK, BQSR, STAT and MERGE 
    Reseqall	Run the whole pipeline of variant calling and phasing, including Varcall, SVcall and Phasing
    Clean		Delete temporary files

    Basic		Run customized steps in `Basicall`
    Reseq		Run customized steps in `Reseqall`
    
#### Config
    python LRTK-SEQ.py Config [options]
    
    Basic options:
    -o --outputdir, output directory.
    
    Advanced options:
    -s --softwarepath, software directory [default: ./bin]
    -d --datasetpath, dataset directory [default: ./dataset]
    
`Config` would generate two configuration files: 1. Basic.config for the components in `Basicall` 2. Reseq.config for the components in `Reseq.config` <br>

#### Basicall
    python LRTK-SEQ.py Basicall [options]
    
    Basic options:
    -i --input, input file for fastq information (Three columns:1.Sample ID; 2.Library ID; 3. Path to the fastqs)
    -o --outputdir, the path to output
    -c --config, configuration file [default: ./config/Basic.config]

    Advanced options:
    -m --minlen, the minimum length (bp) of molecule to be considered [default: 500]
    -p --parallel, the number CPUs can be used in parallel [default: 1] 
    -N --noBX, generate additional fastqs without barcode tags [default: yes]

`Basicall` includes five submodules: `ALN`, `MARK`, `BQSR`, `STAT` and `MERGE`, which would be carried out step by step automatically. Users can also choose any steps to run independently by using `‘Basic’`. <br><br> output file?

eg. <br>
python LRTK.py Basicall -i fqlist.txt -o out_dir -c ./config/Basic.config <br>

##### Basic:
    python LRTK-SEQ.py Basic \<command\> [options]
    
    ALN   generate clean fastq files and correct barcode errors
    MARK      barcode-aware PCR duplicates removal
    BQSR      base quality score recalibration
    STAT      calculate QC statistics, including Cf, Cr, MuFL, NFP etc
    MERGE     merge all bam files from the same sample

###### ALN
    python LRTK-SEQ.py Basic ALN
    
    Basic options:
    -i --input, input file for fastq information (Three columns:1.Sample ID; 2.Library ID; 3. Path to the fastqs)
    -o --outputdir, the path to output
    -c --config, configuration file [default: ./config/Basic.config]

    Advanced options
    -p --parallel, the number CPUs can be used in parallel [default: 1] 
    -N --noBX, generate additional fastqs without barcode tags [default: yes]

`CFQ_ALN` is similar as what `basic` and `align` in [`Long ranger`](https://support.10xgenomics.com/genome-exome/software/pipelines/latest/what-is-long-ranger) do, including barcode error correction, barcode white-listing, attaching barcodes to reads, and read alignment by `bwa`. The description of output files are written in out_dir/Result_list/Basic_ALN_result.txt. Meanwhile, it also generates the input file for the next step `MARK`: out_dir/Result_list/Basic_MARK_input.txt. <br><br>

eg. <br>
python LRTK-SEQ.py Basic ALN -i fqlist.txt -o out_dir -c ./config/Basic.config <br>
###### MARK
    python LRTK.py Basic MARK
    
    Basic options:
    -i --input, input file for BAM information (Three columns: 1.Sample Id; 2.Library Id; 3.Path to BAM, generated by `CFQ_ALN`)
    -o --outputdir, the path to output
    -c --config, configuration file [default: ./config/Basic.config]
    
    Advanced option:
    -p --parallele, the number CPUs can be used in parallel [default: 1] 


`MARK` would merge all the bams from the same library and perform barcode-aware PCR duplicates removal using [`picard`]. The description of output files are written in out_dir/Result_list/Basic_MARK_result.txt. Meanwhile, it also generates the input file for the next step `BQSR`: out_dir/Result_list/Basic_BQSR_input.txt. <br><br>

e.g. <br>
python LRTK-SEQ.py Basic MARK -i out_dir/Result_list/Basic_MARK_input.txt -o out_dir -c ./config/Basic.config <br>



###### BQSR
    python LRTK-SEQ.py Basic BQSR
    
    Basic options:
    -i --input, input file for BAM information (Three columns: 1.Sample Id; 2.Library Id; 3.Path to BAM, generated by `MARK`)
    -o --outputdir, the path to output
    -c --config, configuration file [default: ./config/Basic.config]
    
    Advanced option:
    -p --parallele, the number CPUs can be used in parallel [default: 1] 

`BQSR` performs base quality score recalibration in [`GATK`](https://software.broadinstitute.org/gatk/).

The description of output files are written in out_dir/Result_list/Basic_BQSR_result.txt. Meanwhile, it also generates the input file for the next step `STAT`: out_dir/Result_list/Basic_STAT_input.txt. <br><br>

e.g.<br>
python LRTK-SEQ.py Basic BQSR -i out_dir/Result_list/Basic_BQSR_input.txt -o outputdir -c ./config/Basic.config <br>


###### STAT
    python LRTK-SEQ.py Basic STAT
    
    Basic options
    -i --input, input file for BAM information (Three columns: 1.Sample Id; 2.Library Id; 3.Path to BAM, generated by `MARK`)
    -o --outputdir, the path to output
    -c --config, configuration file [default: ./config/Basic.config]

    
    Advanced options
    -p --parallele, the number CPUs can be used in parallel [default: 1] 
    -m --minlen, the minimum length (bp) of molecule to be considered [default: 500]

The short reads align to the adjacent locations and share the same barcodes should come from the same long DNA fragment. `LRTK-SEQ` reconstructs these long DNA fragments and writes their information into `fragment.csv`. It contains nine columns: 1. fragment id 2. barcode sequence 3. chromosome 4. start position 5. end position 6. fragment length 7. number of reads 8. total length of all the reads in this fragment 9. depth coverage. 

`LRTK-SEQ` also generates four histograms: 1. Unweighted fragment length distribution 2. Weighted fragment length distribution 3. number of fragments per droplet 4. The distribution of sequencing depth per fragment.

The description of output files are written in out_dir/Result_list/Basic_STAT_result.txt. Meanwhile, it also generates the input file for the next step `MERGE`: out_dir/Result_list/Basic_MERGE_input.txt. <br><br>

e.g. <br>
python LRTK-SEQ.py Basic STAT -i out_dir/Result_list/Basic_STAT_input.txt -o out_dir -c ./config/Basic.config <br>



###### MERGE
    python LRTK-SEQ.py Basic MERGE
    
    Basic options:
    -i --input, input file for BAM information (Three columns: 1.Sample Id; 2.Library Id; 3.Path to BAM, generated by `MARK`)
    -o --outputdir, the path to output
    -c --config, configuration file [default: ./config/Basic.config]
    
    Advanced option:
    -p --parallele, the number CPUs can be used in parallel [default: 1] 
    
The description of output files are written in out_dir/Result_list/Basic_MERGE_result.txt. Meanwhile, it also generates the input file for the next module `Reseq`: out_dir/Result_list/Reseq_Varcall_Phasing_input.txt. <br><br>


e.g. <br>
python LRTK-SEQ.py Basic MERGE -i out_dir/Result_list/Basic_MERGE_input.txt -o outputdir -c ./config/Basic.config <br>

#### Reseqall
    python LRTK-SEQ.py Reseqall [options]
    
    Basic options:
    -i --input, the BAM information files from Basicall or ALN
    -o --outputdir, the path to output
    -c --config, configuration file [default: ./config/Reseq.config]
 
    Advanced optionL
    -p --parallel, the number CPUs can be used in parallel [default: 1]
    
 `Basicall` includes five submodules: `ALN`, `MARK`, `BQSR`, `STAT` and `MERGE`, which would be carried out step by step automatically. Users can also choose any steps to run independently by using `‘Basic’`. <br><br> output file?


`Reseqall` includes three submodules: `Varcall`, `SVcall` and `Phasing`, which would be carried out step by step automatically. Users can also choose any steps to run independently by using `‘Reseq’`. output?

eg. <br>
python LRTK-SEQ.py Reseqall -i out_dir/Result_list/Reseq_Varcall_Phasing_input.txt -o out_dir -c ./config/Reseq.config <br>

##### Reseq
    python LRTK-SEQ.py Reseq [options]
    
    Varcall     call SNVs and Indels by GATK
    SVcall      call structure variantions by NAIBR 
    Phasing     phasing variants by HapCUT2

###### Varcall
    python LRTK-SEQ.py Reseq Varcall
    
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
    python LRTK-SEQ.py Reseq SVcall
    
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
    python LRTK-SEQ.py Reseq Phasing
    
    Basic options:
    -i --input, the BAM information files from Basicall or ALN
    -v --vcf, unphased vcf file generated by Varcall
    -o --outputdir, the path to output
    -c --config, configuration file [default: ./config/Reseq.config]

`Phasing` phases the variants in vcf files by [HapCUT2](https://github.com/vibansal/HapCUT2). The phased vcf is written in out_dir/XX.vcf.<br><br>

e.g. <br>
python LRTK-SEQ.py Reseq Phasing -i out_dir/Result_list/Reseq_Varcall_Phasing_input.txt -v out_dir/Result_list/Reseq_Phasing_input.txt -o out_dir -c ./config/Reseq.config <br>

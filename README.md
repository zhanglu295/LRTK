# LRTK
This is the genomic resequencing pipeline for read cloud data.
## Requirements
LRTK was written by python 3, and required packages are: os, sys, gzip, getopt, time, re, subprocess, random, string, glob, collections, multiprocessing, pysam

## Basic usage
python LRTK.py \<command\> [options]

### Command:
Config		Generate configuration file
Basicall		Execute the whole pipeline of basic data preparation
Reseqall		Execute the whole pipeline of resequencing
Denovoall	Execute the whole pipeline of de novo assembly
Clean		delete temporary files

Basic		Execute selected steps for basic data preparation
Reseq		Execute selected steps for resequencing
Denovo		Execute selected steps for de novo assembly

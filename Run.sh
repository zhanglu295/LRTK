python /home/meijunpu/LRTK/LRTK-master/run-LRTK.py QC CFQ -i /home/meijunpu/work_daliy/test_LRTK/fq.txt -o $PWD/result -p 3 1>log.o 2>log.e
python /home/meijunpu/LRTK/LRTK-master/run-LRTK.py QC ALN -i /home/meijunpu/work_daliy/test_LRTK/result/clean_fq.txt -o $PWD/result -p 3 1>log.o 2>log.e
python /home/meijunpu/LRTK/LRTK-master/run-LRTK.py QC MAK -i /home/meijunpu/work_daliy/test_LRTK/result/original_bam.txt -o /home/meijunpu/work_daliy/test_LRTK/result -p 3 1>log.o 2>log.e
python /home/meijunpu/LRTK/LRTK-master/run-LRTK.py QC STAT -i /home/meijunpu/work_daliy/test_LRTK/result/merge_marked_bam.txt -o $PWD/result -p 3


python /home/meijunpu/LRTK/LRTK-master/run-LRTK.py QCall -i /home/meijunpu/work_daliy/test_LRTK/fq.txt -o $PWD/TZ -p 3


python /home/meijunpu/LRTK/LRTK-master/run-LRTK.py Reseq Varcall -i /home/meijunpu/work_daliy/test_LRTK/result/merge_marked_bam.txt -o /home/meijunpu/work_daliy/test_LRTK/result -L /home/meijunpu/work_daliy/test_LRTK/chrList.txt -p 3
python /home/meijunpu/LRTK/LRTK-master/run-LRTK.py Reseq Phasing -i /home/meijunpu/work_daliy/test_LRTK/result/merge_marked_bam.txt -v /home/meijunpu/work_daliy/test_LRTK/result/unphased.vcf.txt -o /home/meijunpu/work_daliy/test_LRTK/result -p 3


python /home/meijunpu/LRTK/LRTK-master/run-LRTK.py Reseqall -i /home/meijunpu/work_daliy/test_LRTK/result/merge_marked_bam.txt -o /home/meijunpu/work_daliy/test_LRTK/result -L /home/meijunpu/work_daliy/test_LRTK/chrList.txt -p 3

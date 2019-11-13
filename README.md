# code for scTAPS

scTAPS
  Call methylation CpG sites/CH sites by Methyl_taps

  *The input bam file must by index.
>python Methy_taps/main.py -h

>Inspecting and pre-processing inputs
>>usage: main.py [-h] 

                -i INPUT_BAM -o OUTPUT_DIR [-p PREFIX] -r REFERENCE
               [-l {directional,non-directional}] [-c {All,CpG,CHH,CHG}]
               [-m {mCtoT,CtoT}] [-k] [-bq MINBQ] [-mq MINMQ]
               [--ignore_orphans] [--ignore_wrong_insert] [--Dip_C]
               [--clip_5_end CLIP_5_END] [--clip_3_end CLIP_3_END]
               [--read_len READ_LEN] [--compress COMPRESS]
               [--N_threads THREADS]

>>The code to call the modifications in the bam/sam reads. Update in 2019.06.11
by Bailey.

>optional arguments:

>>  -h, --help            show this help message and exit

>>  -i INPUT_BAM, --input_bam INPUT_BAM
                        input bam format file containing sequencing reads.
                        (default: None)
                        
>>  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        The output directory to save files. (default: None)
                        
>>  -p PREFIX, --prefix PREFIX
                        The prefix such as sample_id. (default: None)
                        
>>  -r REFERENCE, --reference REFERENCE
                        Reference DNA sequence in FASTA format, such as
                        genome.fa (default: None)
                        
>>  -l {directional,non-directional}, --library {directional,non-directional}
                        The type of library preparation (Default non-
                        directional). (default: non-directional)
                        
>>  -c {All,CpG,CHH,CHG}, --context {All,CpG,CHH,CHG}
                        Explains which cytosine sequence contexts are to be
                        expected in the output file. Default behaviour is all,
                        which includes CpG, CHG, CHH contexts and their sub-
                        contexts for downstream filtering and analysis.
                        (default: All)
                        
>>  -m {mCtoT,CtoT}, --modified_method {mCtoT,CtoT}
                        Specify sequencing method, possible options are CtoT
                        (unmodified cytosines are converted to thymines,
                        bisulfite sequencing-like) and mCtoT (modified
                        cytosines are converted to thymines, TAPS-like).
                        (default: mCtoT)
                        
>>  -k, --skip_clip_overlap
                        Skipping the random removal of overlapping bases
                        between paired-end reads. Not recommended for paired-
                        end libraries, unless the overlaps are removed prior
                        to calling. (Default False) (default: False)
                        
>>  -bq MINBQ, --minimum_base_quality MINBQ
                        Set the minimum base quality for a read base to be
                        used in the pileup (Default 13). (default: 1)
                        
>>  -mq MINMQ, --minimum_mapping_quality MINMQ
                        Set the minimum mapping quality for a read to be used
                        in the pileup (Default 0). (default: 30)
                        
>>  --ignore_orphans      Ignore reads not in proper pairs (Default False).
                        (default: False)
                        
>>  --ignore_wrong_insert
                        Ignore ignore wrong insert size(Default False).
                        (default: False)
                        
>>  --Dip_C               Process reads from 3C genome sequence. (default:
                        False)
                        
>>  --clip_5_end CLIP_5_END
                        remove <int> bp from the 5' end of read. (Default 10)
                        (default: 10)
  
>>  --clip_3_end CLIP_3_END
                        remove <int> bp from the 3' end of read. (Default 10)
                        (default: 10)
  
>>  --read_len READ_LEN   remove <int> bp from the 3' end of read when the
                        fragment length or single read length is less than seq
                        read length (default: 150)
  
>>  --compress COMPRESS, -z COMPRESS
                        Indicates whether the mods file output will be
                        compressed with gzip (Default False). (default: None)
                        
>>  --N_threads THREADS, -t THREADS
                        The number of threads to spawn (Default 1). (default:
                        1)



High-throughtput scTAPS



>sep barcode first

>python sep_scripts/Tn5.barcode.find.8.py -h

>usage: Tn5.barcode.find.8.py [-h] 

                             [-R1 READ1_FASTQ] [-R2 READ2_FASTQ]
                             [-b BARCODE_SEQ_FILE] [-o OUTDIR]
                             [--min_dis MIN_DIS] [--Tn5 TN5]
                             [--Tn5_min_dis TN5_MIN_DIS] [--prefix PREFIX]

The code is used to seperate reads by different cell barcode. Update in
2019.07.01 by Bailey. Copyright (c) Bailey

optional arguments:
  -h, --help            show this help message and exit
  -R1 READ1_FASTQ, --Read1_fastq READ1_FASTQ
                        Read1 fastq file path.
  -R2 READ2_FASTQ, --Read2_fastq READ2_FASTQ
                        Read2 fastq file path.
  -b BARCODE_SEQ_FILE, --barcode_seq_file BARCODE_SEQ_FILE
                        barcode seq ref file: barcode_id barcode_seq.
  -o OUTDIR, --outdir OUTDIR
                        The outdir of output files
  --min_dis MIN_DIS     min dist of (index+barcode seq) and seq (Default 3).
  --Tn5 TN5             The Tn5 seq after barcode(Default
                        AGATGTGTATAAGAGACAG)..
  --Tn5_min_dis TN5_MIN_DIS
                        min dist of Tn5.
  --prefix PREFIX       The out name



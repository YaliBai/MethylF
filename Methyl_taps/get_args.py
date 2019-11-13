import argparse
import os
import sys



def get_args(args=None):
	if args is None:
		args = sys.argv[1:]
	parser = argparse.ArgumentParser(
	description="The code to call the modifications in the bam/sam reads. Update in 2019.05.10 by Bailey.",
	formatter_class = argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-i", "--input_bam",
	required=True,
	help="input bam format file containing sequencing reads.")
	parser.add_argument("-o", "--output_dir",
	type=str,
	required=True,
	help="The output directory to save files.")
	parser.add_argument("-p","--prefix",
	help="The prefix such as sample_id.")
	parser.add_argument("-r", "--reference",
	dest = 'reference',
	type=str,
	required=True,
	help="Reference DNA sequence in FASTA format, such as genome.fa")
	parser.add_argument("-l", '--library',
	choices=['directional', 'non-directional'],
	default='non-directional',
	help="The type of library preparation (Default non-directional).")
	parser.add_argument("-c", "--context",
	choices=['All', 'CpG', 'CHH', 'CHG'],
	default='All',
	help="Explains which cytosine sequence contexts are to be expected in the output file. Default behaviour is all, which includes CpG, CHG, CHH contexts and their sub-contexts for downstream filtering and analysis.")
	parser.add_argument("-m", "--modified_method",
	choices=['mCtoT', 'CtoT'],
	default='mCtoT',
	help="Specify sequencing method, possible options are CtoT (unmodified cytosines are converted to thymines, bisulfite sequencing-like) and mCtoT (modified cytosines are converted to thymines, TAPS-like).")

	parser.add_argument("-k", "--skip_clip_overlap",
	dest = 'skip_clip_overlap',
	action="store_true",
	help="Skipping the random removal of overlapping bases between paired-end reads. Not recommended for paired-end libraries, unless the overlaps are removed prior to calling. (Default False)")
	parser.add_argument('-bq', "--minimum_base_quality",
	dest = 'minBQ',
	default=1,
	help="Set the minimum base quality for a read base to be used in the pileup (Default 13).")
	parser.add_argument('-mq', "--minimum_mapping_quality",
	dest = 'minMQ',
	default =30,
	help="Set the minimum mapping quality for a read to be used in the pileup (Default 0).")

	parser.add_argument('--ignore_orphans',
	dest = 'ignore_orphans',
	action="store_true",
	help="Ignore reads not in proper pairs (Default False).")
	
	parser.add_argument('--ignore_wrong_insert',
	dest = 'ignore_wrong_insert',
        action="store_true",
        help="Ignore ignore wrong insert size(Default False).")


	parser.add_argument('--Dip_C',
	action="store_true",
	help="Process reads from 3C genome sequence.")

	parser.add_argument('--clip_5_end',
	type=int,
	default=10,
	help="remove <int> bp from the 5' end of read. (Default 10)")

	parser.add_argument('--clip_3_end',
	type=int,
	default=10,
	help="remove <int> bp from the 3' end of read.  (Default 10)")

	parser.add_argument('--read_len',
        type=int,
        default=150,
        help="remove <int> bp from the 3' end of read when the fragment length or single read length is less than seq read length")

	parser.add_argument('--compress', '-z',
	help="Indicates whether the mods file output will be compressed with gzip (Default False).")
	parser.add_argument('--N_threads', '-t',
	dest = 'threads',
	type=int,
	default =1,
	help="The number of threads to spawn (Default 1).")




	print('\nInspecting and pre-processing inputs')
	args = vars(parser.parse_args(args))

	if(args['output_dir'][-1] == '/'):
		args['output_dir'] = args['output_dir'][0:-1]
	if not os.path.exists(args['output_dir']):
		os.makedirs(args['output_dir'])
#	if not os.path.exists(args['output_dir'] + '/plots'):
#		os.makedirs(args['output_dir'] + '/plots')

	check_pipeline_input(args)

	print('Done.')
	#return vars(parser.parse_args(args))
	return args


def check_pipeline_input(args):
	assert os.path.exists(args['input_bam']), \
		'Cannot find bam file %s' % args['input_bam']
	assert os.path.exists(args['reference']), \
		'Cannot find reference file %s' % args['reference']
	if not args['prefix']:
		args['prefix']=os.path.basename(os.path.splitext(args['input_bam'])[0])

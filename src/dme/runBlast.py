from dme.Blast import Blast
from dme.settings import *
import argparse

cpus = os.cpu_count()

def main(args):
	blast_obj = Blast(args.input_sequence,
		output_file=args.output_file,
		program=args.program,
		num_threads=args.threads,
		outfmt=0,
		local_database=args.local_database
	)
	blast_obj.run()

def create_parser():
	parser = argparse.ArgumentParser(prog="dme blast", description="{} - {} - BLAST".format(APP_NAME, SOFTWARE_VERSION))
	parser.add_argument('-i','--input_sequence', dest="input_sequence", required=True, help='input file must be in either FASTA (contig and protein) or gzip format! e.g myFile.fasta, myFasta.fasta.gz')
	parser.add_argument('-o','--output_file', dest="output_file", required=True, help="output folder and base filename")
	parser.add_argument('-n','--num_threads', dest="threads", type=int, default=cpus, help="number of threads (CPUs) to use in the BLAST search (default={})".format(cpus))
	parser.add_argument('--program', dest="program", type=str.lower, default="blastp", \
	choices=['blastp','blastn', 'blastx', 'tblastn','tblastx'], required=False, help="blast program to use default: blastp")
	parser.add_argument('--local', dest="local_database", action="store_true", help="use local database (default: uses database in executable directory)")
	parser.add_argument('--debug', dest="debug", action="store_true", help="debug mode")
	return parser

def run():
	parser = create_parser()
	args = parser.parse_args()
	main(args)

if __name__ == '__main__':
	run()

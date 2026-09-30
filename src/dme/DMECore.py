import subprocess
import glob
import time
import gzip
import bz2
import filetype

from Bio import SeqIO

from dme.DMEBase import DMEBase
from dme.Database import Database
from dme.Blast import Blast
from dme.Diamond import Diamond
from dme.ORF import ORF
from dme.Filter import Filter
from dme.Output import Output
from dme.settings import *


class DMECore(DMEBase):
    """Class to predict drug metabolizing gene from protein or nucleotide data based on HMDM detection models."""

    def __init__(self, input_type='contig', input_sequence=None, threads=32, output_file=None, loose=False,
                 clean=True, data='na', aligner='blast', galaxy=None, local_database=False, low_quality=False, debug=False, split_prodigal_jobs=False, include_nudge=False):
        """Creates DMECore object for prediction."""

        o_f_path, o_f_name = os.path.split(os.path.abspath(output_file))

        self.input_type = input_type.lower()
        self.input_sequence = os.path.abspath(input_sequence)
        self.threads = threads
        self.num_sequences = 1
        self.output_file = os.path.abspath(output_file)
        self.loose = loose
        self.clean = clean
        self.data = data
        self.aligner = aligner.lower()
        self.database = galaxy
        self.low_quality = low_quality

        self.local_database = local_database
        self.db = path
        self.dp = data_path

        if self.local_database:
            self.db = LOCAL_DATABASE
            self.dp = LOCAL_DATABASE

        self.working_directory = o_f_path
        self.blast_results_xml_file = ''
        self.debug = debug
        self.split_prodigal_jobs = split_prodigal_jobs
        self.include_nudge = include_nudge
        self.umcompressed_file = ""

        if self.debug:
            logger.setLevel(10)

        super(DMEBase, self).__init__()

    def __repr__(self):
        """Returns DMECore class full object."""
        return f"DMECore({self.__dict__}"

    @classmethod
    def from_string(cls, cmd_string):
        """Creates DMECore object from string."""
        input_type, input_sequence, threads, num_sequences, output_file, aligner, database = cmd_string.split(
            '@')
        return cls(input_type, input_sequence, threads, num_sequences, output_file, aligner, database)

    @classmethod
    def from_args(cls, *initial_data, **kwargs):
        """Creates DMECore object from args."""
        for dictionary in initial_data:
            for key in dictionary:
                if key in ['input_type', 'loose', 'clean', 'aligner']:
                    setattr(cls, key, dictionary[key].lower())
                setattr(cls, key, dictionary[key])

        for key in kwargs:
            if key in ['input_type', 'loose', 'clean', 'aligner']:
                setattr(cls, key, kwargs[key].lower())
            setattr(cls, key, kwargs[key])

        return cls()

    def validate_inputs(self):
        """Validate inputs.

                - validate input file name and out file name
                - validation for mutually exclusive options e.g. protein sequence for contig input_type etc
        """
        if not os.path.exists(self.input_sequence):
            logger.error("input file does not exist: %s", self.input_sequence)
            sys.exit()

        # otherwise you blow up your input when deleting intermediate files
        if self.output_file == self.input_sequence and self.clean:
            logger.error(
                "output path same as input, must specify different path when cleaning to prevent accidental deletion of input files")
            sys.exit()

        kind = filetype.guess(self.input_sequence)
        # logger.info(f"{self.input_sequence} => {kind}")

        if kind is None:
            if self.is_fasta() == False:
                logger.error("invalid fasta")
                sys.exit()
        else:
            if kind.extension in ["gz", "bz2"]:
                if self.is_fasta(kind.extension) == False:
                    logger.error("invalid fasta")
                    sys.exit()
                # uncompressed input and use uncompressed file
                filename = os.path.basename(self.input_sequence)
                umcompressed_file = os.path.join(
                    self.working_directory, f"{filename}.temp.uncompressed.fsa")
                with open(umcompressed_file, "w") as file_out:
                    if kind.extension == "gz":
                        with gzip.open(self.input_sequence, "rt") as handle:
                            file_out.write(handle.read())
                    else:
                        with bz2.open(self.input_sequence, "rt") as handle:
                            file_out.write(handle.read())

                self.input_sequence = umcompressed_file
                self.umcompressed_file = umcompressed_file
            else:
                logger.error(
                    "Sorry, no support for file format %s", kind.mime)
                sys.exit()

        if self.threads > os.cpu_count():
            logger.error("Argument num_threads illegal value, expected (>=1 and =<%s):  given `%s`)",
                         os.cpu_count(), self.threads)
            sys.exit()

    def is_fasta(self, extension=""):
        """Checks for valid fasta format."""
        if extension == "":
            with open(self.input_sequence, "r") as handle:
                fasta = SeqIO.parse(handle, "fasta")
                self.check_record(fasta)
                return True
        elif extension in ["gz", "bz2"]:
            if extension == "gz":
                with gzip.open(self.input_sequence, "rt") as handle:
                    fasta = SeqIO.parse(handle, "fasta")
                    self.check_record(fasta)
            else:
                with bz2.open(self.input_sequence, "rt") as handle:
                    fasta = SeqIO.parse(handle, "fasta")
                    self.check_record(fasta)
            return True
        else:
            return False

    def check_record(self, fasta):
        # check each record in the file
        for record in fasta:
            if any(record.id) == False or any(record.seq) == False:
                return False
            if self.input_type == "contig":
                return self.is_dna(record.seq)
            if self.input_type == "protein":
                return self.is_protein(record.seq)

    # TODO: test
    def _is_fasta(self):
        """Checks for valid fasta format using seqkit stats"""
        cmd = f"seqkit stats --tabular {self.input_sequence} --out-file - | grep -v format"
        result = subprocess.check_output(cmd, shell=True)
        result_dict = result.strip().decode().split("\t")

        if self.input_type == "contig":
            if result_dict[1] == "FASTA" and result_dict[2] == ["DNA"]:
                return True
            else:
                return False
        if self.input_type == "protein":
            if result_dict[1] == "FASTA" and result_dict[2] == ["Protein"]:
                return True
            else:
                return False

    @staticmethod
    def is_dna(sequence):
        #  dna codes
        nucleotide_dict = {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'N': 0, 'U': 0,
                           #  other dna codes
                           'W': 0,  # W = A or T
                           'S': 0,  # S = C or G
                           'M': 0,  # M = A or C
                           'K': 0,  # K = G or T
                           'R': 0,  # R = A or G
                           'Y': 0,  # Y = C or T
                           'B': 0,  # B = C, G, or T
                           'D': 0,  # D = A, G, or T
                           'H': 0,  # H = A, C, or T
                           'V': 0  # V = A, C, or G
                           }

        for base in sequence:
            try:
                nucleotide_dict[base.upper()] += 1
            except Exception as e:
                logger.error("invalid nucleotide fasta due to: %s", e)
                return False
        # logger.info(f"valid nucleotide fasta: '{json.dumps(nucleotide_dict)}'")
        return True

    @staticmethod
    def is_protein(sequence):
        amino_acids_dict = {
            # common symbols between protein and dna codes
            'A': 0, 'T': 0, 'G': 0, 'C': 0, 'N': 0, 'U': 0,
            # other amino acids
            'R': 0, 'D': 0, 'Q': 0, 'E': 0, 'G': 0, 'H': 0, 'I': 0,
            'L': 0, 'K': 0, 'M': 0, 'F': 0, 'P': 0, 'S': 0, 'T': 0,
            'W': 0, 'Y': 0, 'V': 0, 'X': 0, 'Z': 0, 'J': 0, 'B': 0
        }
        count = 0
        for amino_acid in sequence:
            try:
                amino_acids_dict[amino_acid.upper()] += 1
            except Exception as e:
                logger.error("invalid protein fasta due to: %s", e)
                return False

        for a in amino_acids_dict.keys():
            if a not in 'ATGCNU':
                count = count + amino_acids_dict[a]

        if count == 0:
            logger.error("invalid protein fasta: %s", amino_acids_dict)
            return False

        logger.info("valid protein fasta: %s", amino_acids_dict)
        return True

    def __set_xml_filepath(self, fp):
        """Sets blast xml filepath."""
        self.blast_results_xml_file = fp

    def create_databases(self):
        """Creates databases."""
        db_obj = Database(self.local_database)
        db_obj.build_databases()

    def run(self):
        """Runs DMECore."""
        t0 = time.time()
        self.validate_inputs()
        self.create_databases()
        self.run_blast()
        self.filter_process()
        # logger.info("Output......")
        # self.out()
        lapse_time = round(time.time() - t0, 3)
        logger.info('Total running time %s s', lapse_time)

    def clean_files(self):
        """Cleans temporary files."""
        if self.clean == True:
            basename_output_file = os.path.splitext(
                os.path.basename(self.output_file))[0]
            logger.info("Cleaning up temporary files...%s",
                        basename_output_file)
            # remove uncompressed input file
            if self.umcompressed_file != "":
                self.remove_file(self.umcompressed_file)
            # clean working_directory
            self.clean_directory(self.working_directory, basename_output_file)
            d_name, f_name = os.path.split(self.output_file)
            # clean destination_directory
            self.clean_directory(d_name, basename_output_file)
        else:
            logger.info("Clean up skipped.")

    def clean_directory(self, directory, basename_output_file):
        """Cleans files in directory."""
        logger.info(directory)
        files = glob.glob(os.path.join(directory, "*"))
        for f in files:
            if os.path.basename(self.input_sequence) + ".temp" in f and os.path.isfile(f):
                self.remove_file(f)
            if os.path.basename(self.input_sequence) + ".fai" in f and os.path.isfile(f):
                self.remove_file(f)
            if "temp.kmers.xml" in f and os.path.isfile(f):
                self.remove_file(f)
            # if os.path.basename(f)[:3] == "tmp" in f and os.path.isfile(f) and ".temp." in f:
            # self.remove_file(f)
            # if ".temp.directory" in f and os.path.isdir(f):
            # logger.info("Removed directory: {}".format(f))
            # shutil.rmtree(f)

    def remove_file(self, f):
        """Removes file."""
        if os.path.exists(f):
            try:
                logger.info("Removed file: %s", f)
                os.remove(f)
            except Exception as e:
                raise e
        else:
            logger.warning("Missing file: %s", f)

    def out(self):
        """Writes tab-delimited, ggf3 output files."""
        tab_obj = Output(self.output_file)
        tab_obj.run()

    def run_blast(self):
        """Runs blast."""
        if self.input_type == "protein":
            self.process_protein()
        elif self.input_type == "contig":
            self.process_contig()
        else:
            logger.error("Invalid input_type: %s ", self.input_type)
            sys.exit()

    def set_xml_filepath(self, fp):
        """Sets blast xml filepath."""
        # logger.info("set blast xml file: [%s]", fp)
        self.blast_results_xml_file = fp

    def process_protein(self):
        """Process protein sequence(s)."""
        file_name = os.path.basename(self.input_sequence)
        xml_file = os.path.join(self.working_directory,
                                f"{file_name}.temp.blastRes.xml")

        if self.aligner == "diamond":
            diamond_obj = Diamond(self.input_sequence, xml_file,
                                  local_database=self.local_database, num_threads=self.threads)
            diamond_obj.run()
        else:
            blast_obj = Blast(self.input_sequence, xml_file,
                              local_database=self.local_database, num_threads=self.threads)
            blast_obj.run()

        self.set_xml_filepath(xml_file)

    def process_contig(self):
        """Process nuclotide sequence(s)."""
        file_name = os.path.basename(self.input_sequence)
        orf_obj = ORF(input_file=self.input_sequence, threads=self.threads, clean=self.clean,
                      working_directory=self.working_directory, low_quality=self.low_quality, split_prodigal_jobs=self.split_prodigal_jobs)
        orf_obj.contig_to_orf()
        contig_fsa_file = os.path.join(
            self.working_directory, f"{file_name}.temp.contig.fsa")
        blast_results_xml_file = os.path.join(
            self.working_directory, f"{file_name}.temp.contig.fsa.blastRes.xml")

        try:
            if os.stat(contig_fsa_file).st_size > 0:
                # logger.info("work with file %s", contig_fsa_file)
                if self.aligner == "diamond":
                    diamond_obj = Diamond(
                        contig_fsa_file, local_database=self.local_database, num_threads=self.threads)
                    diamond_obj.run()
                else:
                    blast_obj = Blast(
                        contig_fsa_file, local_database=self.local_database, num_threads=self.threads)
                    blast_obj.run()
                self.set_xml_filepath(blast_results_xml_file)
            else:
                self.write_stub_output_file()
                logger.error("no open reading frames (orfs) found.")
        except Exception as e:
            self.write_stub_output_file()
            logger.exception("failed to write orf file, %s", e)
        else:
            # logger.info("success procession orf file")
            pass

    def write_stub_output_file(self):
        # write empty output file if there are no open reading frames
        with open(os.path.join(self.output_file), 'w') as fout:
            fout.write(json.dumps({}))

    # @profile
    def filter_process(self):
        # logger.info("run filter")
        """Filter each detection models and predict resistome(s)."""
        filter_obj = Filter(self.input_type, self.loose, self.input_sequence, self.blast_results_xml_file,
                            os.path.join(self.dp, "hmdm.json"), os.path.basename(self.input_sequence), self.output_file, self.threads, self)
        filter_obj.run()

    def output(self): pass

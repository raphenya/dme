import subprocess
from dme.settings import *
from pathlib import Path


class Database(object):
    """Class to create BLAST databases from a hmdm.json file."""

    def __init__(self, local_database=False):
        """Creates Database object."""
        self.local_database = local_database
        self.db = path
        self.data = data_path
        self.stdout = "2>&1 >> /dev/null"  # "2> /dev/null"

        if self.local_database:
            self.db = LOCAL_DATABASE
            self.data = LOCAL_DATABASE

    def __repr__(self):
        """Returns Database class full object."""
        return "Database({}".format(self.__dict__)

    def build_databases(self):
        """Build BLAST and DIAMOND databases."""
        self.write_fasta_from_json()
        self.write_fasta_from_json("dna")
        # exit("debug...")
        self.make_blast_database_protein("protein")
        self.make_blast_database_dna("dna")
        self.make_blast_database_strains("strains")
        self.make_diamond_database("protein")

    def make_blast_database_protein(self, dbtype="protein"):
        """Build BLAST database from a FASTA file."""
        if os.path.isfile(os.path.join(self.db, "{}db.fsa".format(dbtype))) == True and os.path.exists(os.path.join(self.db, "{}db.fsa".format(dbtype))) == True  \
           and os.path.exists(os.path.join(self.db, "{}.db.phr".format(dbtype))) == True and os.path.exists(os.path.join(self.db, "{}.db.pin".format(dbtype))) == True \
           and os.path.exists(os.path.join(self.db, "{}.db.psq".format(dbtype))) == True:
            logger.info(f"blast DB ({dbtype}) exists.")
            pass
        else:
            type = "prot"
            logger.info(f"create blast DB ({dbtype}).")

            infasta = os.path.join(self.db, f"{dbtype}db.fsa")
            outres = os.path.join(self.db, f"{dbtype}.db")
            os.system(
                f'makeblastdb -in {infasta} -dbtype {type} -out {outres} {self.stdout}')
            # response = subprocess.run(
            #     ["makeblastdb", "-in", infasta, "-dbtype", type, "-out", outres, self.stdout])

    def make_blast_database_dna(self, dbtype="dna"):
        """Build BLAST database from a FASTA file."""
        if os.path.isfile(os.path.join(self.db, "{}db.fsa".format(dbtype))) == True and os.path.exists(os.path.join(self.db, "{}db.fsa".format(dbtype))) == True  \
           and os.path.exists(os.path.join(self.db, "{}.db.nhr".format(dbtype))) == True and os.path.exists(os.path.join(self.db, "{}.db.nin".format(dbtype))) == True \
           and os.path.exists(os.path.join(self.db, "{}.db.nsq".format(dbtype))) == True:
            logger.info(f"blast DB ({dbtype}) exists.")
            pass
        else:
            type = "nucl"
            logger.info(f"create blast DB ({dbtype}).")

            infasta = os.path.join(self.db, f"{dbtype}db.fsa")
            outres = os.path.join(self.db, f"{dbtype}.db")
            os.system(
                f'makeblastdb -in {infasta} -dbtype {type} -out {outres} {self.stdout}')
            # response = subprocess.run(
            #     ["makeblastdb", "-in", infasta, "-dbtype", type, "-out", outres, self.stdout])

    def make_blast_database_strains(self, dbtype="strains"):
        """Build BLAST database from a FASTA file."""

        large_db = False

        # get file size for strains fasta
        if os.path.isfile(os.path.join(self.db, "{}db.fsa".format(dbtype))):
            file_path = Path(os.path.join(self.db, "{}db.fsa".format(dbtype)))
            # Get size in bytes and convert to GiB (binary GB)
            size_bytes = file_path.stat().st_size
            size_gb = size_bytes / (1024 ** 3)

            if size_gb >= 1.0:
                large_db = True

            logger.info(
                f"strainsdb is size_gb: {size_gb:.2f}")

        if os.path.isfile(os.path.join(self.db, "{}db.fsa".format(dbtype))) == True and os.path.exists(
                os.path.join(self.db, "{}db.fsa".format(dbtype))) == True \
                and os.path.exists(os.path.join(self.db, "{}.db.nto".format(dbtype))) == True and os.path.exists(
            os.path.join(self.db, "{}.db.not".format(dbtype))) == True \
                and os.path.exists(os.path.join(self.db, "{}.db.ntf".format(dbtype))) == True \
                and os.path.exists(os.path.join(self.db, "{}.db.njs".format(dbtype))) == True \
                and os.path.exists(os.path.join(self.db, "{}.db.ndb".format(dbtype))) == True \
        and os.path.exists(os.path.join(self.db, "{}.db.nal".format(dbtype))) == True and large_db == True:
            logger.info(f"blast large_db({large_db}) DB ({dbtype}) exists.")
            pass
        elif os.path.isfile(os.path.join(self.db, "{}db.fsa".format(dbtype))) == True and os.path.exists(os.path.join(self.db, "{}db.fsa".format(dbtype))) == True  \
                and os.path.exists(os.path.join(self.db, "{}.db.nhr".format(dbtype))) == True and os.path.exists(os.path.join(self.db, "{}.db.nin".format(dbtype))) == True \
                and os.path.exists(os.path.join(self.db, "{}.db.nsq".format(dbtype))) == True and large_db == False:
            logger.info(f"blast large_db({large_db}) DB ({dbtype}) exists.")
            pass
        else:
            logger.info(f"create blast large_db({large_db}) DB ({dbtype}).")
            type = "nucl"
            infasta = os.path.join(self.db, f"{dbtype}db.fsa")
            outres = os.path.join(self.db, f"{dbtype}.db")
            os.system(
                f'makeblastdb -in {infasta} -dbtype {type} -out {outres} {self.stdout}')

    def make_diamond_database(self, dbtype="protein"):
        """Build DIAMOND database from a FASTA file."""
        if os.path.isfile(os.path.join(self.db, "{}db.fsa".format(dbtype))) == True \
                and os.path.exists(os.path.join(self.db, "{}db.fsa".format(dbtype))) == True and os.path.exists(os.path.join(self.db, "{}.db.dmnd".format(dbtype))) == True:
            logger.info(f"diamond DB ({dbtype}) exists.")
            pass
        else:
            logger.info(f"create diamond DB ({dbtype}).")
            infasta = os.path.join(self.db, f"{dbtype}db.fsa")
            outres = os.path.join(self.db, f"{dbtype}.db")
            os.system(
                f'diamond makedb --quiet --in {infasta} --db {outres} {self.stdout}')
            # response = subprocess.run(["diamond makedb","--quiet","--in",infasta, "--db", outres, self.stdout])

    def make_custom_db(self, in_file, out_file, db_type="nucl", program="blast"):
        if program == 'blast':
            os.system('makeblastdb -in {in_file}  -dbtype {db_type} -out {out_file} \
				{stdout}'.format(in_file=in_file, db_type=db_type, out_file=out_file, stdout=self.stdout))
        else:
            exit("Only NCBI BLAST is supported.")

    def write_fasta_from_json(self, dbtype="protein"):
        """Creates a fasta file from hmdm.json file."""
        if os.path.isfile(os.path.join(self.db, "{}db.fsa".format(dbtype))):
            # logger.info("Database already exists.")
            return
        else:
            try:
                with open(os.path.join(self.data, "hmdm.json"), 'r') as jfile:
                    j = json.load(jfile)
            except Exception as e:
                logger.error(e)
                exit()

            with open(os.path.join(self.db, "{}db.fsa".format(dbtype)), 'w') as fout:
                for i in j:
                    if i.isdigit():
                        # model_type: protein homolog model
                        if j[i]['model_type_id'] == 54:
                            try:
                                pass_bit_score = j[i]['model_param']['blastp_bit_score']['param_value']
                            except KeyError:
                                logger.warning("No bitscore for model (%s, %s). DME will omit this model and keep running."
                                               % (j[i]['model_id'], j[i]['model_name']))
                                logger.info(
                                    "Please let the HMDM Admins know! Email: raphenar@mcmaster.ca")
                            else:
                                try:
                                    for seq in j[i]['model_sequences']['sequence']:
                                        fout.write('>%s_%s | model_type_id: 54 | pass_bitscore: %s | %s\n' % (
                                            i, seq, pass_bit_score, j[i]['HMDM_name']))
                                        fout.write('%s\n' % (
                                            j[i]['model_sequences']['sequence'][seq]['{}_sequence'.format(dbtype)]['sequence']))
                                except Exception as e:
                                    logger.warning("No model sequences for model (%s, %s). DME will omit this model and keep running."
                                                   % (j[i]['model_id'], j[i]['model_name']))
                                    logger.info(
                                        "Please let the HMDM Admins know! Email: raphenar@mcmaster.ca")

    def write_fasta_from_json_dna(self):
        """Creates a fasta file from hmdm.json file."""
        if os.path.isfile(os.path.join(self.db, "dnadb.fsa")):
            # logger.info("Database already exists.")
            return
        else:
            try:
                with open(os.path.join(self.data, "hmdm.json"), 'r') as jfile:
                    j = json.load(jfile)
            except Exception as e:
                logger.error(e)
                exit()

            with open(os.path.join(self.db, "dnadb.fsa"), 'w') as fout:
                for i in j:
                    if i.isdigit():
                        # model_type: protein homolog model
                        if j[i]['model_type_id'] == 54:
                            try:
                                pass_bit_score = j[i]['model_param']['blastp_bit_score']['param_value']
                            except KeyError:
                                logger.warning("No bitscore for model (%s, %s). DME will omit this model and keep running."
                                               % (j[i]['model_id'], j[i]['model_name']))
                                logger.info(
                                    "Please let the HMDM Admins know! Email: raphenar@mcmaster.ca")
                            else:
                                try:
                                    for seq in j[i]['model_sequences']['sequence']:
                                        fout.write('>%s_%s | model_type_id: 54 | pass_bitscore: %s | %s\n' % (
                                            i, seq, pass_bit_score, j[i]['HMDM_name']))
                                        fout.write('%s\n' % (
                                            j[i]['model_sequences']['sequence'][seq]['dna_sequence']['sequence']))
                                except Exception as e:
                                    logger.warning("No model sequences for model (%s, %s). DME will omit this model and keep running."
                                                   % (j[i]['model_id'], j[i]['model_name']))
                                    logger.info(
                                        "Please let the HMDM Admins know! Email: raphenar@mcmaster.ca")

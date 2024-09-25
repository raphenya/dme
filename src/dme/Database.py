from dme.settings import *


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
        self.make_blast_database()
        self.make_blast_database("dna")
        self.make_diamond_database()

    def make_blast_database(self, dbtype="protein"):
        """Build BLAST database from a FASTA file."""
        if os.path.isfile(os.path.join(self.db, "{}db.fsa".format(dbtype))) == True and os.path.exists(os.path.join(self.db, "{}db.fsa".format(dbtype))) == True  \
           and os.path.exists(os.path.join(self.db, "{}.db.phr".format(dbtype))) == True and os.path.exists(os.path.join(self.db, "{}.db.pin".format(dbtype))) == True \
           and os.path.exists(os.path.join(self.db, "{}.db.psq".format(dbtype))) == True:
            logger.info("blast DB exists")
            pass
        else:
            logger.info("create blast DB.")
            type = "prot"
            if dbtype == "dna":
                type = "nucl"
            os.system('makeblastdb -in {} -dbtype {type} -out {} {stdout}'.format(os.path.join(self.db,
                      "{}db.fsa".format(dbtype)), os.path.join(self.db, "{}.db".format(dbtype)), type=type, stdout=self.stdout))

    def make_diamond_database(self, dbtype="protein"):
        """Build DIAMOND database from a FASTA file."""
        if os.path.isfile(os.path.join(self.db, "{}db.fsa".format(dbtype))) == True and os.path.exists(os.path.join(self.db, "{}db.fsa".format(dbtype))) == True \
                and os.path.exists(os.path.join(self.db, "{}.db.dmnd".format(dbtype))) == True:
            logger.info("diamond DB exists")
            pass
        else:
            logger.info("create diamond DB.")
            os.system('diamond makedb --quiet --in {} --db {} {stdout}'.format(os.path.join(self.db,
                      "{}db.fsa".format(dbtype)), os.path.join(self.db, "{}.db".format(dbtype)), stdout=self.stdout))

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

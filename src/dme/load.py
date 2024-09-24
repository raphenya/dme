import shutil
import argparse
from dme.settings import *

'''
loaded_databases example:

{
    "hmdm_canonical": {
      "hmdm_version": "1.0.0"
    },
    "hmdm_variants": {
      "hmdm_version": "1.0.0"
    },
    "hmdm_kmers": {
      "kmer_sizes": ["31","61"]
    }
  }
'''
loaded_databases = {}

# this script is used to load new hmdm.json file to system wide package or local


def get_hmdm_json_version(hmdm_json_path):
    data_version = ""
    if os.path.isfile(hmdm_json_path) == True:
        with open(hmdm_json_path) as json_file:
            json_data = json.load(json_file)
            for item in json_data.keys():
                if item == "version":
                    data_version = json_data[item]
    return data_version


def validate_file(filename):
    try:
        with open(filename) as f:
            out = json.load(f)
            return out
    except ValueError as e:
        logger.error("Invalid json: {}".format(e))
        return None  # or: raise


def get_hmdm_annotation(args_hmdm_annotation, args_local_database):
    hmdm_annotation = args_hmdm_annotation
    # check if hmdm annation is already loaded
    if args_hmdm_annotation is None:
        db = get_location(args_local_database)
        if os.path.isfile(os.path.join(db, "hmdm_reference.fasta")):
            logger.info("hmdm annotation exists")
            hmdm_annotation = os.path.join(db, "hmdm_reference.fasta")
    return hmdm_annotation


def main(args):
    # versions used
    if os.path.isfile(os.path.join(get_location(args.local_database), "loaded_databases.json")):
        # file exists load
        with open(os.path.join(get_location(args.local_database), "loaded_databases.json"), 'r') as fout:
            loaded_databases = json.load(fout)
    else:
        # initialize bwt dict for versions used
        loaded_databases = {
            "hmdm_canonical": {
                "data_version": "N/A"
            },
            "hmdm_variants": {
                "data_version": "N/A"
            },
            "hmdm_kmers": {
                "kmer_sizes": []
            }
        }

    # print args
    if args.debug:
        logger.setLevel(10)
    logger.info(json.dumps(args.__dict__, indent=2))

    if args.hmdm_json is not None:
        # validate json
        if validate_file(args.hmdm_json) == False:
            logger.error("failed to read json file: {}".format(args.hmdm_json))
            exit()
        load_file(args.local_database, args.hmdm_json, "hmdm.json")
        loaded_databases["hmdm_canonical"]["data_version"] = get_hmdm_json_version(
            os.path.join(get_location(args.local_database), "hmdm.json"))

    hmdm_annotation = args.hmdm_annotation

    # if args.hmdm_annotation is not None and (args.wildhmdm_index is None or args.wildhmdm_annotation is None):
    if args.hmdm_annotation is not None:
        load_reference_hmdm_only(
            args.local_database, args.hmdm_annotation, "hmdm_reference.fasta")
    else:
        # check if hmdm annation is already loaded
        hmdm_annotation = get_hmdm_annotation(
            args.hmdm_annotation, args.local_database)

    # and args.hmdm_annotation is not None:
    if args.wildhmdm_index is not None and args.wildhmdm_annotation is not None:
        # load index
        load_file(args.local_database, args.wildhmdm_index,
                  "index-for-model-sequences.txt")
        # load annotation files (hmdm and wildhmdm)
        load_reference_hmdm_and_wildhmdm(
            args.local_database, hmdm_annotation, args.wildhmdm_annotation, "hmdm_wildhmdm_reference.fasta")
        loaded_databases["hmdm_variants"]["hmdm_version"] = args.wildhmdm_version

    if args.kmer_database is not None:
        if args.kmer_size is not None:
            load_file(args.local_database, args.kmer_database,
                      "{}mer_database.json".format(str(args.kmer_size)))
            if args.kmer_size not in loaded_databases["hmdm_kmers"]["kmer_sizes"]:
                loaded_databases["hmdm_kmers"]["kmer_sizes"].append(
                    args.kmer_size)
        else:
            logger.error("Need to specify kmer size when loading kmer files.")

    if args.amr_kmers is not None:
        if args.kmer_size is not None:
            load_file(args.local_database, args.amr_kmers,
                      "amr_{}mer.txt".format(str(args.kmer_size)))
            if args.kmer_size not in loaded_databases["hmdm_kmers"]["kmer_sizes"]:
                loaded_databases["hmdm_kmers"]["kmer_sizes"].append(
                    args.kmer_size)
        else:
            logger.error("Need to specify kmer size when loading kmer files.")

    # and args.hmdm_annotation is not None:
    if args.baits_index is not None and args.baits_annotation is not None:
        logger.info("adding index and fasta for baits")
        # load index
        load_file(args.local_database, args.baits_index,
                  "baits-probes-with-sequence-info.txt")
        # load annotation files (baits)
        load_reference_hmdm_and_baits(
            args.local_database, hmdm_annotation, args.baits_annotation, "hmdm_baits_reference.fasta")

    if args.wildhmdm_index is not None and args.wildhmdm_annotation is not None and args.baits_index is not None and args.baits_annotation is not None:
        # load annotations files for HMDM, VARIANTS and BAITS
        load_reference_hmdm_and_wilhmdm_and_baits(
            args.local_database, hmdm_annotation, args.baits_annotation, args.wildhmdm_annotation, "hmdm_wildhmdm_baits_reference.fasta")
        loaded_databases["hmdm_variants"]["hmdm_version"] = args.wildhmdm_version

    if args.strains_annotation is not None:
        load_file(args.local_database,
                  args.strains_annotation, "strainsdb.fsa")

    # write versions used
    print("write loaded_databases.json to %s",
          os.path.join(get_location(args.local_database)))
    with open(os.path.join(get_location(args.local_database), "loaded_databases.json"), 'w') as fout:
        json.dump(loaded_databases, fout)
    # print out loaded databases
    logger.info(json.dumps(loaded_databases, indent=2))


def load_reference_hmdm_only(local_db, fasta_file, filename):
    load_file(local_db, fasta_file, "hmdm_reference.fasta")
    logger.info("loaded hmdm only for 'dme bwt'.")


def load_reference_baits_only(local_db, fasta_file):
    load_file(local_db, fasta_file, "baits_reference.fasta")
    logger.info("loaded baits only for 'dme bwt'.")


def load_reference_hmdm_and_wilhmdm_and_baits(local_db, hmdm_fasta_file, baits_fasta_file, wildhmdm_fasta_file, filename):
    db = get_location(local_db)
    filenames = []
    filenames.append(hmdm_fasta_file)
    filenames.append(wildhmdm_fasta_file)
    filenames.append(baits_fasta_file)

    # combine the three fastas
    import fileinput
    with open(os.path.join(db, filename), 'w') as fout, fileinput.input(filenames) as fin:
        for line in fin:
            fout.write(line)
    logger.info("loaded hmdm, wildhmdm and baits annotations for 'dme bwt'.")


def load_reference_hmdm_and_baits(local_db, hmdm_fasta_file, baits_fasta_file, filename):
    db = get_location(local_db)

    filenames = []
    filenames.append(hmdm_fasta_file)
    filenames.append(baits_fasta_file)

    # combine the two fastas
    import fileinput
    with open(os.path.join(db, filename), 'w') as fout, fileinput.input(filenames) as fin:
        for line in fin:
            fout.write(line)
    # load baits only file
    load_reference_baits_only(local_db, baits_fasta_file)
    logger.info("loaded hmdm and baits annotations for 'dme bwt'.")


def load_reference_hmdm_and_wildhmdm(local_db, hmdm_fasta_file, wildhmdm_fasta_file, filename):

    db = get_location(local_db)
    # check files
    if hmdm_fasta_file is None:
        # if os.path.isfile(os.path.abspath(hmdm_fasta_file)) == False:
        logger.error("missing hmdm reference file")
        exit()
    if os.path.isfile(os.path.abspath(wildhmdm_fasta_file)) == False:
        logger.error("missing wildhmdm reference file")
        exit()

    filenames = []
    filenames.append(hmdm_fasta_file)
    filenames.append(wildhmdm_fasta_file)

    # combine the two fastas
    import fileinput
    with open(os.path.join(db, filename), 'w') as fout, fileinput.input(filenames) as fin:
        for line in fin:
            fout.write(line)
    logger.info("loaded hmdm and wildhmdm annotations for 'dme bwt'.")


def get_location(local_db):
    db = ""

    # path to save hmdm.json file and database files
    if local_db == True:
        db = LOCAL_DATABASE
        # create directory if it doesn't exist
        if not os.path.exists(LOCAL_DATABASE):
            os.makedirs(LOCAL_DATABASE)
    else:
        db = data_path

    return db


def load_file(local_db, filepath, filename, validate_json=False):

    db = get_location(local_db)

    try:
        # copy new file
        if filename in ["strainsdb.fsa"] and local_db == False:
            # copy file to _db directory and not _data directory
            shutil.copyfile(filepath, os.path.join(path, filename))
        else:
            shutil.copyfile(filepath, os.path.join(db, filename))
            logger.info("file {} loaded ok".format(filename))
    except Exception as e:
        logger.warning("failed to copy json file: {}".format(e))


def create_parser():
    parser = argparse.ArgumentParser(
        prog="dme load", description="{} - {} - Load".format(APP_NAME, SOFTWARE_VERSION))
    parser.add_argument('-i', '--hmdm_json', required=False,
                        help='must be a hmdm database json file')

    parser.add_argument('--hmdm_annotation', required=False,
                        help="annotated reference FASTA")

    parser.add_argument('--wildhmdm_annotation',
                        required=False, help="annotated reference FASTA")
    parser.add_argument('--wildhmdm_index', required=False,
                        help="wildhmdm index file (index-for-model-sequences.txt)")
    parser.add_argument('--wildhmdm_version', required=False,
                        help="specify variants version used")

    parser.add_argument('--baits_annotation', required=False,
                        help="annotated reference FASTA")
    parser.add_argument('--baits_index', required=False,
                        help="baits index file (baits-probes-with-sequence-info.txt)")

    parser.add_argument('--kmer_database', required=False,
                        help="json of kmer database")
    parser.add_argument('--amr_kmers', required=False,
                        help="txt file of all amr kmers")
    parser.add_argument('--kmer_size', required=False,
                        help="kmer size if loading kmer files")

    parser.add_argument('--strains_annotation',
                        required=False, help="annotated strain k-mers")

    parser.add_argument('--local', dest="local_database", action="store_true",
                        help="use local database (default: uses database in executable directory)")
    parser.add_argument('--debug', dest="debug",
                        action="store_true", help="debug mode")
    return parser


def run():
    parser = create_parser()
    args = parser.parse_args()
    main(args)


if __name__ == "__main__":
    run()

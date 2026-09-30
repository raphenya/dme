import shutil
import argparse
from dme.settings import *
import tempfile
from bs4 import BeautifulSoup as bs
import requests
import json
import re


def main(args):
    debug = ""
    if args.debug:
        logger.setLevel(10)
        debug = "--debug"

    local_database = ""
    if args.local_database:
        local_database = "--local"

    hmdm_cannonical_version, hmdm_variants_version = get_versions()

    logger.info(json.dumps(args.__dict__, indent=2))
    logger.info("hmdm cannonical version: {}".format(hmdm_cannonical_version))
    logger.info("hmdm variants version: {}".format(hmdm_variants_version))

    # Create the directory
    directory = tempfile.mkdtemp(
        prefix=os.path.join(os.getcwd(), "dme_autoload_"))
    print("Directory '%s' created" % directory)
    print("=================================== DOWNLOAD HMDM CANONICAL DATA ===================================")
    # get latest hmdm database
    data = os.path.join(directory, "data")
    hmdm_data = os.path.join(directory, "hmdm_data")
    os.system("wget -O {data} --no-check-certificate https://hmdm.mcmaster.ca/download/0/data-v{hmdm_cannonical_version}.tar.bz2".format(
        data=data,
        hmdm_cannonical_version=hmdm_cannonical_version
    )
    )
    os.system("mkdir -p {hmdm_data}".format(hmdm_data=hmdm_data))
    os.system(
        "tar xf {data} -C {hmdm_data}".format(data=data, hmdm_data=hmdm_data))

    print("=================================== DOWNLOAD HMDM VARIANTS DATA ===================================")
    variants = os.path.join(directory, "variants")
    hmdm_variants = os.path.join(directory, "hmdm_variants")
    os.system("wget -O {variants} --no-check-certificate https://hmdm.mcmaster.ca/download/6/prevalence-v{hmdm_variants_version}.tar.bz2".format(
        variants=variants,
        hmdm_variants_version=hmdm_variants_version
    )
    )
    os.system("mkdir -p {hmdm_variants}".format(hmdm_variants=hmdm_variants))
    os.system("tar xf {variants} -C {hmdm_variants}".format(
        variants=variants, hmdm_variants=hmdm_variants))
    os.system(
        "gunzip {hmdm_variants}/*.gz".format(hmdm_variants=hmdm_variants))

    print("=================================== HMDM CANONICAL ANNOTATIONS ===================================")
    os.system(
        "dme hmdm_annotation --input {hmdm_data}/hmdm.json".format(hmdm_data=hmdm_data))

    print("=================================== HMDM VARIANTS ANNOTATIONS ===================================")
    os.system("dme hmdm_annotation --input_directory {hmdm_variants} --version {hmdm_variants_version} --hmdm_json {hmdm_data}/hmdm.json".format(
        hmdm_variants=hmdm_variants,
        hmdm_variants_version=hmdm_variants_version,
        hmdm_data=hmdm_data
    )
    )

    print("=================================== CLEAN OLD DATABASES ===================================")
    os.system("dme clean {debug} {local_database}".format(
        local_database=local_database, debug=debug))

    print("=================================== LOAD DATABASES ===================================")
    os.system("dme load \
	--hmdm_json {hmdm_data}/hmdm.json \
	--hmdm_annotation hmdm_database_v{hmdm_cannonical_version}.fasta \
	--wildhmdm_index {hmdm_variants}/index-for-model-sequences.txt \
	--wildhmdm_version {hmdm_variants_version} \
	--wildhmdm_annotation wildhmdm_database_v{hmdm_variants_version}.fasta \
	--kmer_database {hmdm_variants}/61_kmer_db.json \
	--amr_kmers {hmdm_variants}/all_amr_61mers.txt \
	--kmer_size 61 \
	{local_database} {debug}".format(
        hmdm_data=hmdm_data,
        hmdm_cannonical_version=hmdm_cannonical_version,
        hmdm_variants=hmdm_variants,
        hmdm_variants_version=hmdm_variants_version,
        local_database=local_database,
        debug=debug
    )
    )

    print("=================================== CHECK LOADED DATABASES ===================================")
    os.system(
        "dme database -v --all {local_database}".format(local_database=local_database))

    if args.clean:
        print("=================================== CLEAN UP ===================================")
        os.system("rm {data}".format(data=data))
        os.system("rm {variants}".format(variants=variants))
        os.system("rm {hmdm_data}/* ".format(hmdm_data=hmdm_data))
        os.system("rm {hmdm_variants}/* ".format(hmdm_variants=hmdm_variants))
        os.system("rm -r {hmdm_data} ".format(hmdm_data=hmdm_data))
        os.system("rm -r {hmdm_variants}".format(hmdm_variants=hmdm_variants))
        os.system("rm -r {}".format(directory))
        os.system("rm hmdm_database_v{hmdm_cannonical_version}.fasta".format(
            hmdm_cannonical_version=hmdm_cannonical_version))
        os.system("rm wildhmdm_database_v{hmdm_variants_version}.fasta".format(
            hmdm_variants_version=hmdm_variants_version))

    print("=================================== DONE ===================================")


def get_versions():
    r = requests.get('https://hmdm.mcmaster.ca/download')
    soup = bs(r.content, 'lxml')
    data = [item['href'] if item.get('href') is not None else item['src']
            for item in soup.select('[href^="/download/0"]')]
    data_version = valid_version(
        re.search(r'\s*([\d.].([\d.]).([\d.]))', data[0]).group(1))
    prev = [item['href'] if item.get('href') is not None else item['src']
            for item in soup.select('[href^="/download/6"]')]
    prev_version = valid_version(
        re.search(r'\s*([\d.].([\d.]).([\d.]))', prev[0]).group(1))
    return (data_version, prev_version)


def valid_version(s):
    msg = "Not a version: '{0}'.".format(s)
    pattern = re.compile('^\s*([\d.]).([\d.]).([\d.])$')
    if pattern.match(s) is not None:
        return s
    else:
        raise argparse.ArgumentTypeError(msg)


def create_parser():
    parser = argparse.ArgumentParser(
        prog="dme auto_load", description="{} - {} - Automatic Load".format(APP_NAME, SOFTWARE_VERSION))
    parser.add_argument('--local', dest="local_database", action="store_true",
                        help="use local database (default: uses database in executable directory)")
    parser.add_argument('--clean', dest="clean",
                        action="store_true", help="removes temporary files")
    parser.add_argument('--debug', dest="debug",
                        action="store_true", help="debug mode")
    return parser


def run():
    parser = create_parser()
    args = parser.parse_args()
    main(args)


if __name__ == "__main__":
    run()

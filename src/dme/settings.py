import os
import sys
import logging
import json

from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio import SeqIO

# ====================================================================================
# FUNTIONS
# ====================================================================================


def determine_path():
    try:
        root = __file__
        if os.path.islink(root):
            root = os.path.realpath(root)
        return os.path.dirname(os.path.abspath(root))
    except:
        print("I'm sorry, but something is wrong.")
        print("There is no __file__ variable. Please contact the author.")
        sys.exit()

# ====================================================================================
# FILEPATHS
# ====================================================================================


script_path = determine_path()

path = os.getenv('DB_PATH', os.path.join(script_path, "_db/"))
data_path = os.getenv('DATA_PATH', os.path.join(script_path, "_data/"))

# ====================================================================================
# LOGGING CONFIG
# ====================================================================================
level = logging.WARNING
logger = logging.getLogger(__name__)
logger.setLevel(level)

# detailed log
formatter = logging.Formatter(
    '%(levelname)s %(asctime)s : (%(filename)s::%(funcName)s::%(lineno)d) : %(message)s')
# basic log
# formatter = logging.Formatter('%(levelname)s %(asctime)s : %(message)s')

stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)

logger.addHandler(stream_handler)

LOCAL_DATABASE = os.path.join(os.getcwd(), "localDB")

APP_NAME = "Drug Metabolising Enzyme (DME)"
SOFTWARE_VERSION = "1.0.2"
SOFTWARE_SUMMARY = 'Use the Drug Metabolising Enzyme (DME) to predict drug-metabolizing enzymes from protein or nucleotide \
data based on homology models. Check https://hmdm.mcmaster.ca/download for software and data updates. \
Receive email notification of monthly HMDM updates via the HMDM Mailing List \
(https://mailman.mcmaster.ca/mailman/listinfo/hmdm-l)'

GALAXY_PROJECT_WRAPPER = 'GALAXY_DATABASE must contain the following: \
data files (hmdm.json, proteindb.fsa),\
diamond blast database (protein.db.dmnd), \
ncbi blast database (protein.db.phr, protein.db.pin, protein.db.psq)'

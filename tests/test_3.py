import pytest
import os
import subprocess as sp
from dme.MainBase import MainBase

inputs = "inputs/"
outputs = "outputs/"
working_directory = os.getcwd()

# Run all tests with
# pytest test_3.py -v -rxs --color=auto --durations=0
# or
# pytest test_3.py -v -rxs --color=auto --durations=0 -k "create"


@pytest.fixture
def dme():
    return MainBase(api=True)


def test_create_local_with_annotations_db(dme):
    parser = dme.load_args()
    f = os.path.join(working_directory, inputs, "{}".format("hmdm.json"))
    annotation = os.path.join(
        working_directory, inputs, "{}".format("hmdm_database_v1.0.0.fasta"))
    fsa = os.path.join(working_directory, inputs, "{}".format("strainsdb.fsa"))
    db = os.path.join(working_directory, "localDB", "{}".format("hmdm.json"))
    strains = os.path.join(working_directory, "localDB",
                           "{}".format("strainsdb.fsa"))
    ref = os.path.join(working_directory, "localDB",
                       "{}".format("hmdm_reference.fasta"))
    dme.load_run(parser.parse_args([
        '--hmdm_json', f,
        '--strains_annotation', fsa,
        '--hmdm_annotation', annotation,
        '--local',
        '--debug'
    ]))

    assert (os.path.isfile(f) and os.path.exists(db)) and (
        os.path.isfile(fsa) and os.path.exists(strains)) and (os.path.isfile(annotation) and os.path.exists(ref)) == True


def run_dme(dme, read_one, read_two, aligner, threads, output_file):
    parser = dme.bwt_args()
    dme.bwt_run(parser.parse_args([
        '--read_one', read_one,
        '--read_two', read_two,
        '--aligner', aligner,
        '--threads', threads,
        '--output_file', output_file,
        '--clean',
        '--debug',
        '--local'
    ]))


def check_gene(gene):
    if gene == "Collinsella aerofaciens hypothetical protein (COLAER_01846)":
        return True
    else:
        return False


def test_dme_bwt_kma(dme):
    read_one = os.path.join(working_directory, inputs,
                            "{}".format("R1.fastq"))
    read_two = os.path.join(working_directory, inputs,
                            "{}".format("R2.fastq"))
    filename = "output_bwt_kma"
    output_file = os.path.join(
        working_directory, outputs, "{}".format(filename))

    run_dme(dme, os.path.join(working_directory, inputs, read_one), os.path.join(
        working_directory, inputs, read_two), "kma", "3", output_file)

    gene = sp.getoutput("cat {gene_summary_output} | cut -f1 | grep '{gene}'".format(
        gene_summary_output="{}.gene_mapping_data.txt".format(output_file),
        gene="Collinsella aerofaciens hypothetical protein (COLAER_01846)"
    ))
    assert check_gene(gene) == True

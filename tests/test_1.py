import pytest
import os
import json
import csv
from dme.MainBase import MainBase

inputs = "inputs/"
outputs = "outputs/"
alignment_tool = "blast"
working_directory = os.getcwd()

# Run all tests with
# pytest test_1.py -v -rxs --color=auto --durations=0
# or
# pytest test_1.py -v -rxs --color=auto --durations=0 -k "protein"


@pytest.fixture
def dme():
    return MainBase(api=True)


def run_dme(dme, input_type, input_sequence, output_file):
    parser = dme.main_args()
    dme.main_run(parser.parse_args([
        '--input_type', input_type,
        '--input_sequence', input_sequence,
        '--output_file', output_file,
        '--alignment_tool', alignment_tool,
        '--clean',
        '--include_loose',
        '--include_nudge',
        '--low_quality',
        '--debug',
        '--local'
    ]))


def validate_results(filepath, perc_identity=0, name='', type_match=''):
    pi = ""
    name = ""
    tm = ""
    filename = os.path.basename(filepath)
    f = os.path.join("{}".format(filepath))
    if os.path.isfile(f):
        with open(f) as json_file:
            json_data = json.load(json_file)
            for i in json_data:
                if i not in ["_metadata"]:
                    for j in json_data[i]:
                        for k in json_data[i][j]:
                            pi = json_data[str(i)][str(j)]["perc_identity"]
                            name = json_data[str(i)][str(j)]["HMDM_name"]
                            tm = json_data[str(i)][str(j)]["type_match"]
                        if pi == perc_identity and name == name and tm == type_match:
                            # print(pi, name, tm)
                            return True
            return False
    else:
        print("missing file: {}".format(f))
        return False


def test_dme_protein_sequence(dme):

    filename = "test-prot.fasta"
    output_file = os.path.join(
        working_directory, outputs, f"{filename}.json")
    run_dme(dme, 'protein', os.path.join(
        working_directory, inputs, filename), output_file)
    print("output_file: {}".format(output_file))
    assert validate_results(output_file, 100, 'AcbK', 'Perfect') == True

import pytest
import os
import json
from dme.MainBase import MainBase

inputs = "inputs/"
outputs = "outputs/"
alignment_tool = "blast"
working_directory = os.getcwd()

# Run all tests with
# pytest test_2.py -v -rxs --color=auto --durations=0
# or
# pytest test_2.py -v -rxs --color=auto --durations=0 -k "drug"


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
    drug_name = ""
    tm = ""
    f = os.path.join("{}".format(filepath))
    if os.path.isfile(f):
        with open(f) as json_file:
            json_data = json.load(json_file)
            for i in json_data:
                if i not in ["_metadata"]:
                    for j in json_data[i]:
                        for k in json_data[i][j]:
                            pi = json_data[str(i)][str(j)]["perc_identity"]
                            for _, category in json_data[str(i)][str(j)]["HMDM_category"].items():
                                if category["category_hmdm_name"] == name and category["category_hmdm_class_name"] == "Drug":
                                    drug_name = category["category_hmdm_name"]
                            tm = json_data[str(i)][str(j)]["type_match"]
                        if pi == perc_identity and drug_name == name and tm == type_match:
                            return True
            return False
    else:
        print("missing file: {}".format(f))
        return False


def test_dme_drug_match(dme):

    filename = "test-nucl.fasta"
    output_file = os.path.join(
        working_directory, outputs, f"{filename}.json")
    run_dme(dme, 'contig', os.path.join(
        working_directory, inputs, filename), output_file)

    assert validate_results(output_file, 100, 'Levodopa', 'Perfect') == True

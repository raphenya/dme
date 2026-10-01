import pytest
import os
from dme.MainBase import MainBase

inputs = "inputs/"
outputs = "outputs/"
alignment_tool = "diamond"
working_directory = os.getcwd()

# Run all tests with
# pytest test_0.py -v -rxs --color=auto --durations=0
# or
# pytest test_0.py -v -rxs --color=auto --durations=0 -k "create"


@pytest.fixture
def dme():
    return MainBase(api=True)


def test_create_local_db(dme):
    parser = dme.load_args()
    f = os.path.join(working_directory, inputs, "{}".format("hmdm.json"))
    fsa = os.path.join(working_directory, inputs, "{}".format("strainsdb.fsa"))
    db = os.path.join(working_directory, "localDB", "{}".format("hmdm.json"))
    strains = os.path.join(working_directory, "localDB",
                           "{}".format("strainsdb.fsa"))
    dme.load_run(parser.parse_args([
        '--hmdm_json', f,
        '--strains_annotation', fsa,
        '--local',
        '--debug'
    ]))

    assert (os.path.isfile(f) and os.path.exists(db)) and (
        os.path.isfile(fsa) and os.path.exists(strains)) == True

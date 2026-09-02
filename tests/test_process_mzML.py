import sys
import os
import subprocess
import json
import glob
from deepdiff import DeepDiff

def eprint(*args, **kwargs): print(*args, file=sys.stderr, **kwargs)

from runassessor.mzML_assessor import MzMLAssessor


# Set up to run assess_mzMLs.py on the test data file 'chlo_6_tiny.mzML.gz'
tests_dir = os.path.dirname(os.path.abspath(__file__))
print(tests_dir)
project_root = os.path.dirname(tests_dir)
print(project_root)

src_dir = os.path.join(project_root, "src")
assess_mzMLs_module = "runassessor.runassessor"

tests_data_dir = os.path.join(tests_dir, "data")
input_mzML_file = "chlo_6_tiny.mzML.gz" # Located in tests/data directory which is the current working directory for subprocess
input_mzML_filepath = os.path.join(tests_data_dir, input_mzML_file)

# The current output of 'assess_mzMLs.py' is 'study_metadata.json' in the tests/data directory
# The expected output is 'chlo_6_tiny_study_metadata.json' in the tests/expected_output directory
current_study_metadata = os.path.join(tests_data_dir, "study_metadata.json")
tests_expected_output_dir = os.path.join(tests_dir, "expected_output")
expected_study_metadata = os.path.join(tests_expected_output_dir, "chlo_6_tiny_study_metadata.json")


def test_assess_mzMLs():

    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(filter(None, [src_dir, env.get("PYTHONPATH")]))

    output = subprocess.run(
        [sys.executable, '-m', assess_mzMLs_module, '--verbose', input_mzML_file],
        cwd=tests_data_dir,
        env=env,
        capture_output=True,
        text=True,
        check=True
    )

    # Testing verbosity
    stdout_lines = []
    for line in output.stdout.splitlines():
        if line != '':
            stdout_lines.append(line)
    print(stdout_lines)

    stderr_lines = []
    for line in output.stderr.splitlines():
        if line != '':
            stderr_lines.append(line)
    eprint(stderr_lines)

    assert len(stdout_lines) == 0
    assert len(stderr_lines) >= 11      # 11 "INFO" lines

    # Make sure the current generated 'study_metadata.json' file exactly matches the expected output 'chlo_6_tiny_study_metadata.json'
    with open(current_study_metadata, 'r') as study_metadata, open(expected_study_metadata, 'r') as chlo_6_tiny_study_metadata:
        current = json.load(study_metadata)
        expected = json.load(chlo_6_tiny_study_metadata)
    
    diff = DeepDiff(expected, current)
    assert not diff, \
        f"Current generated study_metadata.json for '{input_mzML_file}' does not match expected output." \
        f"\n\nDifferences found:\n{diff.pretty()}\n"

    #### If these assertions pass, then go ahead and remove the test output files
    study_metadata_files = os.path.join(tests_data_dir, "study_metadata.*")
    for file_path in glob.glob(study_metadata_files):
        if os.path.exists(file_path):
            os.remove(file_path)


def test_mzML_assessor():

    assessor = MzMLAssessor(mzml_file=input_mzML_filepath, verbose=1)

    model_data = assessor.read_header()
    assert model_data['accession'] == 'MS:1001911'
    assert model_data['category'] == 'pureHCD'
    assert model_data['name'] == 'Q Exactive'

    stats = assessor.read_spectra()
    assert stats['n_spectra'] == 100
    assert stats['n_ms1_spectra'] == 92
    assert stats['n_ms2_spectra'] == 8

    ROIs = assessor.assess_lowend_composite_spectra()
    assert ROIs['TMT6_126']['peak']['mode_bin']['n_spectra'] == 2

    assert assessor.metadata['state']['status'] != 'ERROR'



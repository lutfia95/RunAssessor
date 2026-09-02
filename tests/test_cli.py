import importlib.metadata
import os
import subprocess
import sys


tests_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(tests_dir)
src_dir = os.path.join(project_root, "src")


def run_cli(*arguments):
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(filter(None, [src_dir, env.get("PYTHONPATH")]))
    return subprocess.run(
        [sys.executable, "-m", "runassessor.runassessor", *arguments],
        env=env,
        capture_output=True,
        text=True,
        check=True,
    )


def test_version_matches_distribution_metadata():
    result = run_cli("--version")

    assert result.stdout.strip().endswith(importlib.metadata.version("runassessor"))


def test_help_documents_singular_sdrf_option():
    result = run_cli("--help")

    assert "--write_sdrf_file" in result.stdout
    assert "--write_sdrf_files" not in result.stdout

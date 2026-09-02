# RunAssessor
The RunAssessor package provides functionality for assessing a set of input mass spectrometry proteomics mzML files and extracting all information that can
be gleaned from the files prior to any attempt to identify the peptides therein.


## Installation
Assuming that Python 3.x (3.8 or greater) is already installed, simply install the RunAssessor package:
```bash
pip install runassessor
```
Note that `pip` will try to install a command-line script in the customary Python location for such scripts
and `pip` will ususally warn if the directory is not in the current PATH. For example, under Windows if you
install as an unprivileged user, it might appear in C:\Users\me\AppData\Roaming\Python\Python312\Scripts.

## Quick Start

RunAssessor processes mass spectrometry data in the `.mzML` format. Compressed `.mzML.gz` files are also supported.
```bash
runassessor ~/file_path/mass_spec_data.mzML.gz
```

RunAssessor can also accept multiple `.mzML` or '.mzML.gz' files at a time 
```bash
runassessor --write_sdrf_file --verbose *.mzML
```

### Running a quick test

In the RunAssessor code repository there is a small `.mzML.gz` file that can be used to quickly try RunAssessor
```bash
git clone https://github.com/RunAssessor/RunAssessor
cd RunAssesor
cd tests
cd data
runassessor chlo_6_tiny.mzML.gz
```

There is also a set of larger `.mzML.gz` files of different kinds that can be downloaded and used to test RunAssessor in `RunAssessor/tests/large_files/`.
The files in `large_files` will generate a more substantial JSON and summary file compared to `chlo_6_tiny.mzML.gz`, but take substantially longer to
download and processing, depending on your Internet and computer speed.
```bash
cd tests/large_files
python large_test_file_downloader.py
runassessor *.mzML.gz   # Or single file of your choice
```


# Usage of the main RunAssessor programs

## assess_mzMLs.py

The main `runassessor` program  accepts a set of input mzML files along with various parameters and processes
the mzML files accordingly. This is generally the only program that needs to be executed.
If multiple files are specified, they will generally be run in parallel using as
many cores as are available on your computer. 

### RunAssessor Command-line Arguments

| Argument | Description | Default | Example |
|----------|-------------|---------|---------|
| `--metadata_filepath` | Filename and path of the ouput metadata file. This parameter may be used to change the root name and directory path (relative or absolute) of the output file (affecting the JSON file, the TSV file, optional PDF files and SDRF files) | `study_metadata.json` | `--metadata_filepath ~/PXD123456/PXD123456_metadata.json` |
| `--preserve_existing_metadata` | If set, read the existing metadata file, preserve its contents, and build on it (default is to overwrite current file). | Not set | `--preserve_existing_metadata` |
| `--n_threads` | Number of files to read in parallel. | Number of CPU cores | `--n_threads 4` |
| `--write_fragmentation_type_files` | If set, writes a fragmentation type file for each mzML input. | Not set | `--write_fragmentation_type_files` |
| `--verbose` | If set, prints more detailed processing information. Use multiple times to increase verbosity. | Not set | `--verbose` |
| `--version` | Prints the version of the script/tool. | Not set | `--version` |
| `--write_sdrf_file` | If set, writes an SDRF file based on the metadata text file. | Not set | `--write_sdrf_file` |
| `--include_sdrf_provenance` | If set, then all information written to the SDRF file will have provenance tags (this is experimental non-standard SDRF, not for routine use). | Not set | `--include_sdrf_provenance` |
| `--write_pdfs` | If set, generates a PDF of delta time and ppm graphs from precursor stats, and a PDF of neutral loss windows from composite intensities. | Not set | `--write_pdfs` |
| `--write_ions` | If set, writes a TSV file of the ion three-sigma table, which lists all fragment ions found. | Not set | `--write_ions` |
| `--ms_ms_units` | Sets units for MS/MS tolerance recommendations: mz or ppm (defaults to ppm for HR and mz for LR) | ppm for HR and mz for LR | `--ms_ms_units mz` |


## assess_metadata.py

The assess_metadata.py program reads an existing study_metadata.json file and recomputes recommended parameters and potentially outputs to
different formats. This program is generally not routinely used, but can be helpful if, for example, a large set of input files are analyzed
and a study_metadata.json file is written, but the --write_sdrf_file was not provided to assess_mzMLs.py. This tool can take the previous analysis
and write out the SDRF file based on the existing JSON file, without reprocessing all mzML files again.

### Command-line Arguments
| Argument | Description | Default | Example |
|----------|-------------|---------|---------|
| `--metadata_file` | Specify a metadata file if different than the default. | `study_metadata.json` | `--metadata_file ~/file_path/custom_metadata.json` |
| `--write_sdrf_file` | If set, writes an SDRF file based on the metadata file. Can be passed multiple times. | Not set | `--write_sdrf_file` |
| `--verbose` | If set, enables verbose output. Use multiple times to increase verbosity. | Not set | `--verbose` |


## rename_msruns.py

Renames raw files in current directory to replace special characters in filenames with underscores and writes a mapping document of changes.
Some special characters in filename can cause problems in downstream processing tools (for example: spaces, parentheses, ampersands, and more). This program
fixes the filenames in a repeatable way and records how the filenames have been adjusted in an output file. This is not needed for normal RunAssessor operation,
but may be useful for large-scale reprocessing efforts when certain tools have trouble with special characters in filenames.




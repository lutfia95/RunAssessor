#!/usr/bin/env python3

import sys
import os
import argparse
import multiprocessing
import copy
from datetime import datetime
from importlib.metadata import version as distribution_version
import timeit
from matplotlib.backends.backend_pdf import PdfPages
from pypdf import PdfReader, PdfWriter

from runassessor.mzML_assessor import MzMLAssessor
from runassessor.metadata_handler import MetadataHandler
from runassessor.graph_generator import GraphGenerator

def eprint(*args, **kwargs): print(*args, file=sys.stderr, **kwargs)


####################################################################################################
#### Process one mzML
def process_job(job):

    #### Assess the mzML file
    assessor = MzMLAssessor(job['filename'], metadata=job['metadata'], verbose=job['verbose'], frag_units=job['fragmentation_units'])
    assessor.read_header()
    assessor.read_spectra(write_fragmentation_type_file=job['write_fragmentation_type_files'])
    assessor.assess_lowend_composite_spectra()
    assessor.assess_neutral_loss_composite_spectra()

    if assessor.metadata['state']['status'] != 'ERROR' or 'multiple fragmentation types' in assessor.metadata['state']['message']:
        assessor.assess_ROIs()
        
        
    return assessor



####################################################################################################
#### Main function for command-line usage
def main():

    argparser = argparse.ArgumentParser(description='Read one or more mzML files and extract some knowledge from them')
    argparser.add_argument('--metadata_filepath', action='store', help='Filename and path of the ouput metadata file (defaults to study_metadata.json in the current working directory)')
    argparser.add_argument('--preserve_existing_metadata', action='count', default=0, help='If set, read the existing metadata file, preserve its contents, and build on it (default is to overwrite current file)')
    argparser.add_argument('--n_threads', action='store', type=int, help='Set the number of files to process in parallel (defaults to number of cores)')
    argparser.add_argument('--write_fragmentation_type_files', action='count', help='If set, write a fragmentation_type file for each mzML')
    argparser.add_argument('--verbose', action='count', help='If set, print more information about ongoing processing' )
    argparser.add_argument('--version', action='version', version=f'%(prog)s {distribution_version("runassessor")}')
    argparser.add_argument('--write_sdrf_file', action='count', help='If set, then write an SDRF file based on the metadata file')
    argparser.add_argument('--include_sdrf_provenance', action='count', help='If set, then all information written to the SDRF file will have provenance tags (non-standard SDRF)')
    argparser.add_argument('--write_pdfs', action='count', help='If set, then generate a PDF of delta time and ppm graphs from precursor stats')
    argparser.add_argument('--write_ions', action='count', help='If set, then write a TSV file of fragment ion peaks found')
    argparser.add_argument('--ms_ms_units', choices=['mz', 'ppm'], help='Sets units for MS/MS tolerance recommendations: mz or ppm (defaults to ppm for HR and mz for LR)')
    argparser.add_argument('files', type=str, nargs='+', help='Filenames of one or more mzML files to read')
    params = argparser.parse_args()

    #### Set verbose level
    verbose = params.verbose
    if verbose is None:
        verbose = 0
    else:
        verbose = 1
    if verbose >= 1:
        timestamp = str(datetime.now().isoformat())
        eprint(f"INFO: Launching RunAssessor at {timestamp}")
    t0 = timeit.default_timer()

    #### Loop over all the files to ensure that they are really there before starting work
    n_files = 0
    for file in params.files:
        if not os.path.isfile(file):
            eprint(f"ERROR: File '{file}' not found or not a file")
            return
        n_files += 1
    if verbose:
        eprint(f"INFO: Found {n_files} input files to process")

    #### Initialize the metadata handler
    study = MetadataHandler(params.metadata_filepath, verbose=params.verbose, write_ions=params.write_ions, frag_units=params.ms_ms_units)

    #### If the user wants to read and preserve an existing file, then read or create it
    if params.preserve_existing_metadata > 0:
        result = study.read_or_create()
    else:
        result = study.create()
    if result is None or result != 'OK':
        eprint(f"ERROR: Unable to initialize study metadata information. Halting.")
        return

    #### Set up the list of jobs to process
    jobs = []

    #### Loop over all the files and put them in a list of jobs for parallel processing
    for file in params.files:

        #### If this file is already in the metadata store, purge it and generate new
        if file in study.metadata['files']:
            if params.preserve_existing_metadata > 0:
                if verbose >= 1:
                    eprint(f"INFO: Already have results for '{file}'. Skipping..")
                continue
            else:
                study.metadata['files'][file] = None

        job = {
            'filename': file,
            'metadata': copy.deepcopy(study.metadata),
            'verbose': verbose,
            'write_fragmentation_type_files': params.write_fragmentation_type_files,
            'fragmentation_units': params.ms_ms_units
            }
        jobs.append(job)


    #### Compute how many files in parallel to process
    n_cpus = multiprocessing.cpu_count()
    n_simultaneous_jobs = n_cpus
    if n_files < n_simultaneous_jobs:
        n_simultaneous_jobs = n_files
    n_threads = params.n_threads
    if n_threads is None or n_threads == 0:
        n_threads = n_simultaneous_jobs
    elif n_threads < 0:
        if -1 * n_threads > n_cpus - 1:
            eprint(f"ERROR: Parameter n_threads ({n_threads} too low for number of CPUs ({n_cpus}))")
            return
        n_threads = n_cpus + n_threads
        if n_threads < n_simultaneous_jobs:
            n_threads = n_files


    #### Now process the jobs in parallel
    if n_threads == 1:
        eprint(f"INFO: Processing 1 input file", end='', flush=True)
    else:
        eprint(f"INFO: Processing {n_threads} files in parallel with multiprocessing (one file per CPU)", end='', flush=True)
    pool = multiprocessing.Pool(processes=n_threads)
    results = pool.map_async(process_job, jobs)
    pool.close()
    pool.join()
    eprint('')

    #### Coalesce the results
    results = results.get()
    for result in results:
        for filename in result.metadata['files']:
            study.metadata['files'][filename] = result.metadata['files'][filename]
        study.metadata['state'] = result.metadata['state']
        study.metadata['state']['problems'] = result.metadata['problems']

    #### Infer parameters based on the latest data
    study.infer_search_criteria()
    study.generate_sdrf_table(include_provenance=params.include_sdrf_provenance)

    #### Write out our state of mind
    study.store()

    #### Write out SDRF table if parameter given
    if params.write_sdrf_file:
        key_value_file = study.find_txt_file()
        study.read_txt_file(key_value_file)
        sdrf_filename = study.infer_sdrf_filename()
        study.write_sdrf_file(sdrf_filename) 

    if verbose >= 1:
        timestamp = str(datetime.now().isoformat())
        t1 = timeit.default_timer()
        eprint(f"INFO: RunAssessor finished in {t1-t0:.2f} seconds at {timestamp}")
    
    #### Write out graphs is parameter given
    if params.write_pdfs:
        grapher = GraphGenerator(params.metadata_filepath, verbose=params.verbose)
        metadata_name = grapher.buildGraphs()
        eprint("INFO: Delta Graphs generated and stored")

        nl_pdf = metadata_name.replace(".json", ".NLplots.pdf")
        # Write cover page documentation to pdf
        nl_coverpage = os.path.join(os.path.dirname(os.path.abspath(__file__)), "neutral_loss_windows_documentation.pdf")
        writer = PdfWriter()
        cover_reader = PdfReader(nl_coverpage)
        for page in cover_reader.pages:
            writer.add_page(page)
        # Generate neutral loss windows and add to the pdf
        with PdfPages(nl_pdf) as pdf:
            for assessor in results:
                grapher.plot_precursor_loss_composite_spectra(assessor=assessor, pdf=pdf)
        output_reader = PdfReader(nl_pdf)
        for page in output_reader.pages:
            writer.add_page(page)
        with open(nl_pdf, "wb") as f:
            writer.write(f)
        eprint("INFO: Neutral Loss Graphs generated and stored")


if __name__ == "__main__": main()

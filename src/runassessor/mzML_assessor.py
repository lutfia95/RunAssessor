#!/usr/bin/env python3

import sys
import os
import argparse
import os.path
import timeit
import re
import numpy
import gzip
from lxml import etree
from functools import partial
import math

#### Import technical modules and pyteomics
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from pyteomics import mzml
import warnings
import scipy.optimize as so

from runassessor.metadata_handler import MetadataHandler

numpy.seterr(all='ignore')
warnings.filterwarnings("ignore", category=so.OptimizeWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

def eprint(*args, **kwargs): print(*args, file=sys.stderr, **kwargs)


####################################################################################################
#### mzML Assessor class
class MzMLAssessor:


    ####################################################################################################
    #### Constructor
    def __init__(self, mzml_file, metadata=None, verbose=None, frag_units=None):
        self.mzml_file = mzml_file

        #### Store the provided metadata or create a template to work in
        if metadata is None:
            self.metadata = { 'files': { mzml_file: {} }, 'state': { 'status': 'OK' } }
        else:
            self.metadata = metadata
            metadata['files'][mzml_file] = {}

        #### Create a place to store our composite spectra for analysis
        self.composite = {}

        #### Create a storage area for referenceable param groups, which are not really handled by pyteomics
        self.referenceable_param_group_list = {}

        #### Create a dictionary to hold 3sigma values to find tolerance
        self.all_3sigma_values_away = {}

        #### stores precursor data 
        self.precursor_stats = {}

        #### Set verbosity
        if verbose is None:
            verbose = 0
        self.verbose = verbose

        #### Set up type of units for ms/ms tolerances
        # frag_units can only be mz, ppm, or None
        self.fragmentation_units = frag_units 

        #### Creates a list of supported composite types that the code can handle
        self.supported_composite_type_list = ['HR_HCD', 'LR_IT_CID', 'LR_HCD', 'HR_QTOF', 'HR_IT_ETD']

        #### Values used to calculate tolerance reccomendations
        self.ppm_error = 5
        self.mz_error = 0.3

        self.skipped_spectrum_types = {}

        self.composite_spectrum_attributes = {
                'lowend': {
                    'HR_HCD': {
                        'minimum': 100,
                        'maximum': 400,
                        'binsize': 0.0002
                    },
                    'LR_IT_CID': {
                        'minimum': 100,
                        'maximum': 460,
                        'binsize': 0.05
                    },
                    'LR_HCD': {
                        'minimum': 100,
                        'maximum': 460,
                        'binsize': 0.05
                    },
                    'LR_IT_ETD': {
                        'minimum': 100,
                        'maximum': 400,
                        'binsize': 0.05
                    },
                    'HR_IT_ETD':{
                        'minimum': 100,
                        'maximum': 400,
                        'binsize': 0.05
                    },
                    'HR_QTOF':{
                        'minimum': 100,
                        'maximum': 400,
                        'binsize': 0.001
                    }
                },
                'precursor_loss': {
                    'HR_HCD': {
                        'minimum': 0,
                        'maximum': 100,
                        'binsize': 0.001
                    },
                    'LR_IT_CID': {
                        'minimum': 0,
                        'maximum': 100,
                        'binsize': 0.05
                    },
                    'LR_HCD': {
                        'minimum': 0,
                        'maximum': 100,
                        'binsize': 0.05
                    },
                    'LR_IT_ETD': {
                        'minimum': 0,
                        'maximum': 100,
                        'binsize': 0.05
                    },
                    'HR_IT_ETD': {
                        'minimum': 0,
                        'maximum': 100,
                        'binsize': 0.05
                    },
                    'HR_QTOF':{
                        'minimum': 0,
                        'maximum': 100,
                        'binsize': 0.001
                    }
                }
            }


    ####################################################################################################
    #### Get empty stats
    def get_empty_stats(self):
        empty_stats = {
                'n_spectra': 0,
                'n_ms0_spectra': 0,
                'n_ms1_spectra': 0,
                'n_ms2_spectra': 0,
                'n_ms3_spectra': 0,
                'n_length_zero_spectra': 0,
                'n_length_lt10_spectra': 0,
                'n_HR_HCD_spectra': 0,
                'n_LR_HCD_spectra': 0,
                'n_HR_IT_CID_spectra': 0,
                'n_LR_IT_CID_spectra': 0,
                'n_HR_IT_ETD_spectra': 0,
                'n_LR_IT_ETD_spectra': 0,
                'n_HR_EThcD_spectra': 0,
                'n_LR_EThcD_spectra': 0,
                'n_HR_ETciD_spectra': 0,
                'n_LR_ETciD_spectra': 0,
                'n_HR_QTOF_spectra': 0,
                'n_unknown_fragmentation_type_spectra': 0,
                'high_accuracy_precursors': 'unknown',
                'fragmentation_type': 'unknown',
                'fragmentation_tag': 'unknown'
            }
        return empty_stats


    ####################################################################################################
    #### Read spectra
    def read_spectra(self, write_fragmentation_type_file=None):

        #### Set up information
        t0 = timeit.default_timer()
        stats = self.get_empty_stats()
        ms_one_tolerance_dict = {}
        self.metadata['files'][self.mzml_file]['spectra_stats'] = stats

        #### Store the fragmentation types in a list for storing in the fragmentation_type_file
        scan_fragmentation_list = []

        #### Show information
        if self.verbose >= 1:
            eprint(f"\nINFO: Assessing mzML file {self.mzml_file}", end='', flush=True)
            progress_intro = False

        #### If the mzML is gzipped, then open with zlib, else a plain open
        if self.mzml_file.endswith('.gz'):
            infile = gzip.open(self.mzml_file)
        else:
            infile = open(self.mzml_file, 'rb')

        #### Try storing all spectra for multithreaded processing
        #spectra = []

        #### Read spectra from the file
        with mzml.read(infile) as reader:
            #try:
            if True:
                for spectrum in reader:

                    #### Debugging. Print the data structure of the first spectrum
                    #if stats['n_spectra'] == 0:
                    #    auxiliary.print_tree(spectrum)
                    #    print(spectrum)
                    #    return

                    #### Extract the MS level of the spectrum
                    try:
                        ms_level = spectrum['ms level']
                    #### If no MS level is found, see if it can be extracted from a referenceable param group
                    except:
                        ms_level = 0
                        if 'ref' in spectrum:
                            group_id = spectrum['ref']
                            if group_id in self.referenceable_param_group_list and 'ms level' in self.referenceable_param_group_list[group_id]:
                                ms_level = int(self.referenceable_param_group_list[group_id]['ms level'])

                    #### Look for a filter string and parse it
                    if 'filter string' in spectrum['scanList']['scan'][0]:
                        filter_string = spectrum['scanList']['scan'][0]['filter string']
                        self.parse_filter_string(filter_string,stats)
                        if filter_string is None or filter_string == '':
                            #### MSFragger generated mzML can have empty filter strings
                            self.log_event('WARNING','EmptyFilterLine',"Filter line is present but empty. This can happen with MSFragger-written mzML. Not good.")
                            stats[f"n_ms{ms_level}_spectra"] += 1
                            if ms_level > 1:
                                stats['n_unknown_fragmentation_type_spectra'] += 1


                    #### There's only a filter string for Thermo data, so for others, record a subset of information
                    else:
                        stats[f"n_ms{ms_level}_spectra"] += 1
                        if self.metadata['files'][self.mzml_file]['instrument_model']['category'] == 'QTOF':
                            stats['high_accuracy_precursors'] = 'true'
                            stats['fragmentation_type'] = 'HR_QTOF'
                            stats['fragmentation_tag'] = 'HR QTOF'
                            if ms_level > 1:
                                stats['n_HR_QTOF_spectra'] += 1

                    #### If the ms level is greater than 2, fail
                    if ms_level > 4:
                        self.log_event('ERROR','MSnTooHigh',f"MS level is greater than we can handle at '{spectrum['ms level']}'")
                        break

                    #### If the ms level is 2, then examine it for information
                    if ms_level == 2 and 'm/z array' in spectrum:
                        try:
                            precursor_mz = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']
                        except:
                            precursor_mz = 0.0

                        #### Try to get the charge information
                        try:
                            charge_state = int(spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['charge state'])
                        except:
                            charge_state = 'unknown'
                        charge_stat = f"n_charge_{charge_state}_precursors"
                        if charge_stat not in stats:
                            stats[charge_stat] = 0
                        stats[charge_stat] += 1

                        # Get the scan number for the spectrum
                        nativeId = spectrum['id']
                        scan_number = -1
                        match = re.search(r'scan=(\d+)',nativeId)
                        if match:
                            scan_number = int(match.group(1))
                        else:
                            match = re.search(r'index=(\d+)',nativeId)
                            if match:
                                scan_number = int(match.group(1))

                        #### Extract the isolation window information
                        self.record_isolation_window_stats(spectrum, stats, ms_one_tolerance_dict, scan_number)

                        #self.add_spectrum(spectrum,spectrum_type,precursor_mz)
                        peaklist = {
                            'm/z array': spectrum['m/z array'],
                            'intensity array': spectrum['intensity array']
                        }

                        #### Check for zero length and very sparse spectra
                        if len(spectrum['m/z array']) == 0:
                            stats['n_length_zero_spectra'] += 1
                        if len(spectrum['m/z array']) < 10:
                            stats['n_length_lt10_spectra'] += 1

                        #spectra.append([stats,precursor_mz,peaklist])
                        self.add_spectrum([stats,precursor_mz,peaklist])

                        # If the user requested writing of a fragmentation_type_file, record that information
                        if write_fragmentation_type_file is not None:
                            scan_fragmentation_list.append(f"{scan_number}\t{stats['fragmentation_tag']}")

                    #### Update counters and print progress
                    stats['n_spectra'] += 1
                    if self.verbose >= 1:
                        if stats['n_spectra']/500 == int(stats['n_spectra']/500):
                            if not progress_intro:
                                #eprint("INFO: Reading spectra.. ", end='')
                                progress_intro = True
                            #eprint(f"{stats['n_spectra']}.. ", end='', flush=True)
                            eprint(".", end='', flush=True)
            #except:
            else:
                self.log_event('ERROR','MzMLCorrupt',f"Pyteomics threw an error reading mzML file! File may be corrupt. Check file '{self.mzml_file}'")

        infile.close()
        #if self.verbose >= 1: eprint("")
        
        #### Find the tolerance of the MS1 scan
        self.find_ms_one_tolerance(ms_one_tolerance_dict)

        # If the user requested writing of a fragmentation_type_file, write it out
        if write_fragmentation_type_file is not None:
            filename = self.mzml_file
            filename = re.sub(r'\.gz','',filename)
            filename = re.sub(r'\.mzML','',filename)
            filename += '.fragmentation_types.tsv'
            eprint(f"\nINFO: Writing fragmentation type file {filename}", end='', flush=True)
            with open(filename,'w') as outfile:
                if stats['fragmentation_type'] != 'multiple':
                    scan_fragmentation = scan_fragmentation_list[0]
                    scan_fragmentation = re.sub(r'^\d+','*',scan_fragmentation)
                    outfile.write(scan_fragmentation+'\n')
                else:
                    for scan_fragmentation in scan_fragmentation_list:
                        outfile.write(scan_fragmentation+'\n')

        #### If there are no spectra that were saved, then we're done with this file
        if stats['n_ms2_spectra'] == 0:
            self.log_event('ERROR','NoMS2Scans',f"This mzML file has no MS2 scans. Check file '{self.mzml_file}'")

        #### Try to distinguish DDA vs DIA
        if 'isolation_window_last_scans' in stats:
            window_stats = {}
            for window, window_details in stats['isolation_window_last_scans'].items():
                n_scan_deltas = len(stats['isolation_window_last_scans'][window]['scan_deltas'])
                if n_scan_deltas <= 2:
                    for scan_delta, n_scans in stats['isolation_window_last_scans'][window]['scan_deltas'].items():
                        if scan_delta not in window_stats:
                            window_stats[scan_delta] = 0
                        window_stats[scan_delta] += n_scans
            stats['n_different_isolation_windows'] = len(stats['isolation_window_last_scans'])
            stats['isolation_window_periodicities'] = window_stats
            if len(stats['isolation_window_periodicities']) <= 5: # Really should just be 1 instead of huge, but allow for future creativity
                stats['acquisition_type'] = 'DIA'
            else:
                stats['acquisition_type'] = 'DDA'
            if len(stats['isolation_window_last_scans']) > 500:
                stats['isolation_window_last_scans'] = { "ommitted": "too many isolation windows to report" }
            if len(stats['isolation_window_periodicities']) >= 10:
                stats['isolation_window_periodicities'] = { "ommitted": "too many isolation window periodicities to report" }

        #### Try to distinguish DDA vs DIA
        elif 'isolation_window_full_widths' in stats:
            smallest_window = 999999
            largest_window = 0
            for window_size in stats['isolation_window_full_widths']:
                window_size = float(window_size)
                if  window_size < smallest_window:
                    smallest_window = window_size
                if window_size > largest_window:
                    largest_window = window_size
            if largest_window <= 3.0:
                stats['acquisition_type'] = 'DDA'
            elif largest_window >= 15.0:
                stats['acquisition_type'] = 'DIA'
            else:
                stats['acquisition_type'] = 'ambiguous'
        else:
            stats['acquisition_type'] = 'unknown'


        #### Now that we've read anything, try to process the data multithreaded
        #n_threads = 8
        #if self.verbose >= 1: eprint(f"INFO: processing {len(spectra)} spectra multithreaded with {n_threads} threads")
        #pool = ThreadPool(n_threads)
        #results = pool.imap_unordered(self.add_spectrum, spectra)
        #results = pool.map(self.add_spectrum, spectra)
        #pool.close()
        #pool.join()
        #if self.verbose >= 1: eprint(f"INFO: Done.")

        #### Write the composite spectrum to a file
        #with open('zzcomposite_spectrum.pickle','wb') as outfile:
        #    pickle.dump(self.composite,outfile)

        #### Print final timing information
        t1 = timeit.default_timer()
        eprint(f"\nINFO: Read {stats['n_spectra']} spectra from {self.mzml_file} in {t1-t0:.2f} sec ({stats['n_spectra']/(t1-t0):.2f} spectra per sec)", end='', flush=True)

        return stats

    ####################################################################################################
    #### Finds the tolerance of the MS1 scan
    def find_ms_one_tolerance(self, ms_one_tolerance_dict):


        delta_ppm_list = []
        delta_time_list = []

        #ms_tol_dict format: int(mz): {"mz":[], "time":[]}

        for key in ms_one_tolerance_dict:
            if len(ms_one_tolerance_dict[key]["mz"]) > 1:
                self.find_deltas(ms_one_tolerance_dict[key], delta_ppm_list, delta_time_list)

            else:
                pass
        

        try:
            #Fit delta_time_list with histogram
            counts_ppm, bins_ppm = numpy.histogram(delta_ppm_list, bins=60)
            bin_centers_ppm = 0.5 * (bins_ppm[1:] + bins_ppm[:-1])
        except:
            pass

        try:
            soft_raised_fixed_beta_ppm = partial(self.precursor_raised_gaussian, beta=1.5)
            p0_ppm = [max(counts_ppm), numpy.mean(delta_ppm_list), numpy.std(delta_ppm_list), 0.1]
            fit_ppm, _ = curve_fit(soft_raised_fixed_beta_ppm, bin_centers_ppm, counts_ppm, p0=p0_ppm)
        
            #Find chi^2 value for a goodness of fit
            expected_counts = soft_raised_fixed_beta_ppm(bin_centers_ppm, *fit_ppm)
            chi_squared = numpy.sum(((counts_ppm - expected_counts) ** 2) / numpy.where(expected_counts == 0, 1, expected_counts))
            num_bins = len(counts_ppm)
            num_params = 4
            dof = num_bins - num_params
            chi_squared_red_ppm = chi_squared / dof
        except:
            pass
 

        
        try:
            #Fit delta_time_list with histogram
            counts_time, bins_time = numpy.histogram(delta_time_list, bins=100)
            bin_centers_time = 0.5 * (bins_time[1:] + bins_time[:-1])
        except:
            pass

        try:
            max_peak_index = numpy.argmax(counts_time)
            main_peak = bin_centers_time[max_peak_index]

            counts_after_peak = counts_time[max_peak_index + 2:]
            counts_before_peak = counts_time[:max_peak_index - 2]

            mean_counts_before_peaks = numpy.mean(counts_before_peak)
            mean_counts_after_peaks = numpy.mean(counts_after_peak)
            if numpy.isnan(mean_counts_before_peaks):
                mean_counts_before_peaks = 0.0

            if numpy.isnan(mean_counts_after_peaks):
                mean_counts_after_peaks = 0.0

            

            p0 = [main_peak, mean_counts_before_peaks, max(counts_time), mean_counts_after_peaks, 1.0]
            lower_bounds = [-numpy.inf, -numpy.inf, -numpy.inf, -numpy.inf, 0]   # pulse_duration >= 0, tau >= 0
            upper_bounds = [numpy.inf,numpy.inf, numpy.inf, numpy.inf, numpy.inf]

            fit_time, cov = curve_fit(self.time_exp_decay, bin_centers_time, counts_time, p0=p0, bounds=(lower_bounds, upper_bounds))
            expected_counts = self.time_exp_decay(bin_centers_time, *fit_time)

            #Compute 5 bins before and 5 bins after for sharpness
            t_before = fit_time[0] - 5
            t_after = fit_time[0] + 5

            y_before = self.time_exp_decay(numpy.array([t_before]), *fit_time)[0]
            y_after = self.time_exp_decay(numpy.array([t_after]), *fit_time)[0]

            bound_y_before = max(y_before, 1e-8)
            bound_y_after = max(y_after, 1e-8)

            # Compute chi-squared
            #window_mask = (bin_centers_time >= t_before) & (bin_centers_time <= t_after)
            #counts_window = counts_time[window_mask]
            #expected_window = expected_counts[window_mask]

            # Compute chi-squared within the window
            epsilon = 1e-8
            chi_squared = numpy.sum(((counts_time - expected_counts) ** 2) / (counts_time+epsilon))
            dof = len(counts_time) - len(fit_time)
            chi_squared_red_time = chi_squared / dof
        except:
            pass

        ppm_fit_call = "no fit"
        time_fit_call = "no fit"

        try:
            if chi_squared_red_ppm < 1.5:
                ppm_fit_call = "good fit found, okay to use"
            elif chi_squared_red_ppm > 1.5 and chi_squared_red_ppm < 10:
                ppm_fit_call = "uncertain about fit, manually check ppm graph"
            else:
                ppm_fit_call = "probably not a good fit, suggested tolerance may not be accurate, check ppm graph"
                
            if chi_squared_red_time < 10:
                time_fit_call = "good fit, dynamic exclusion window found"
            elif chi_squared_red_time > 10 and chi_squared_red_time < 100:
                time_fit_call = "uncertain about fit, manually check time graph"
            else:
                time_fit_call = "probably not a good fit, suggested tolerance may not be accurate, check time graph"
        except:
            pass



        try:
            self.precursor_stats["precursor tolerance"] = {
                            "ppm_fit call": ppm_fit_call,
                            "histogram_ppm":{"counts": counts_ppm.tolist(), "bin_edges": bins_ppm.tolist(), "bin_centers": bin_centers_ppm.tolist()}}
        
        except:
           self.precursor_stats['precursor tolerance']["good ppm_fit?"] = ppm_fit_call
           self.precursor_stats['precursor tolerance']["histogram_ppm"] = "no delta_ppm values recorded"

        try:
            
            self.precursor_stats['precursor tolerance']['fit_ppm'] = {"lower_three_sigma (ppm)":-fit_ppm[2]*3, "upper_three_sigma (ppm)":fit_ppm[2]*3, "delta_ppm peak": fit_ppm[1], "intensity": fit_ppm[0], "sigma (ppm)": fit_ppm[2], "y_offset": fit_ppm[3], "chi_squared":chi_squared_red_ppm, "recommended precursor tolerance (ppm)": math.ceil(math.sqrt(self.ppm_error**2 + (fit_ppm[2]*3)**2))}
        
        except:
            self.precursor_stats['precursor tolerance']['fit_ppm'] = "no guassian ppm curve fit"

        try:
            self.precursor_stats["dynamic exclusion window"]= {
                                "dynamic_exclusion_time_fit_call": time_fit_call,
                                "histogram_time": {"counts": counts_time.tolist(),"bin_edges": bins_time.tolist(),"bin_centers": bin_centers_time.tolist()}}
        except:
            self.precursor_stats["dynamic exclusion window"]["good dynamic exclusion time fit?"] = ppm_fit_call
            self.precursor_stats["dynamic exclusion window"]["histogram_time"] = "no delta_time values recorded"

        try:
            self.precursor_stats['dynamic exclusion window']["fit_pulse_time"] = {"pulse start": fit_time[0], "inital level":fit_time[1], "peak level":fit_time[2], "final level":fit_time[3], "decay constant":fit_time[4], "peak to preceding level ratio": fit_time[2]/bound_y_before, "peak to following level ratio": fit_time[2]/bound_y_after, "chi_squared": chi_squared_red_time}
            
        except:
            self.precursor_stats["dynamic exclusion window"]["fit_pulse_time"] = 'no dynamic exclusion window fit'


    ####################################################################################################
    #### Finds the different deltas in a ms_one_tolerance_dict key     
    def find_deltas(self, info_dict, delta_ppm_list, delta_time_list):
        ions = info_dict['mz']
        times = info_dict['time']
        p_a = 0
        p_b = p_a + 1
        p_c = 0
        used_index = set()

        max_time = 100 # Maximum time in seconds to search for the next same precursor
        min_time = 0 # Minimum time in seconds to search for the next same precursor
        max_delta_ppm = 10 # Maximum distance to search for the next same precursor

        if (len(times) != len(ions)) or not len(ions):
            eprint('WARNING: Unable to find tolerance for MS1 scan')

        while p_c < len(ions) and p_b < len(ions): #Figure out what end condition is
            delta_time = (times[p_b] - times[p_a]) * 60
            delta_ppm = ((ions[p_b] - ions[p_a]) / ions[p_a]) * 10**6

            if delta_time >= min_time and delta_time <= max_time and abs(delta_ppm) <= max_delta_ppm:
                delta_ppm_list.append(delta_ppm)
                delta_time_list.append(delta_time)
                used_index.add(p_b)
                p_a = p_b
                p_b = p_a + 1

            else:
                p_b += 1
                if p_b >= len(ions) or (times[p_b] - times[p_a]) * 60 > max_time:
                    used_index.add(p_c)
                    while p_c in used_index:
                         p_c += 1
                    p_a = p_c
                    p_b = p_a + 1


    ####################################################################################################
    #### Returns a gaussian fit for delta_ppm
    def precursor_raised_gaussian(self, x, a, x0, sigma, baseline, beta=1.5):
        delta = numpy.abs(x - x0)
        delta = numpy.clip(delta, 1e-12, None)
        return a * numpy.exp(-delta**beta / (2 * sigma**beta)) + baseline
    
    
    ####################################################################################################
    #### Pulse exp decay function with peak_level param
    def time_exp_decay(self, t, pulse_start, initial_level, peak_level, new_level, tau):
        exponent_rise = numpy.clip(-(t - pulse_start) * 50, -700, 700)
        rise = 1 / (1 + numpy.exp(exponent_rise))

        # peak at pulse_start, no explicit "plateau"
        y = initial_level + (peak_level - initial_level) * rise

        safe_tau = max(tau, 1e-8)
        decay_mask = t >= pulse_start
        decay = numpy.zeros_like(t)
        decay[decay_mask] = (peak_level - new_level) * numpy.exp(-(t[decay_mask] - pulse_start) / safe_tau)
        y[decay_mask] = new_level + decay[decay_mask]

        return y



        

    ####################################################################################################
    #### Record the isolation window information and precursor ion and scan time information
    def record_isolation_window_stats(self, spectrum, stats, ms_one_tolerance_dict, scan_number):
        isolation_window_target_mz = None
        isolation_window_lower_offset = None
        isolation_window_upper_offset = None
        precursor_ion = None
        start_scan_time = None

        if 'isolation_window_full_widths' not in stats:
            stats['isolation_window_full_widths'] = {}
        if 'isolation_window_last_scans' not in stats:
            stats['isolation_window_last_scans'] = {}

        try:
            isolation_window_target_mz = spectrum['precursorList']['precursor'][0]['isolationWindow']['isolation window target m/z']
            isolation_window_lower_offset = spectrum['precursorList']['precursor'][0]['isolationWindow']['isolation window lower offset']
            isolation_window_upper_offset = spectrum['precursorList']['precursor'][0]['isolationWindow']['isolation window upper offset']
        except Exception as e:
            return

        if isolation_window_lower_offset is not None and isolation_window_upper_offset is not None:
            full_width = isolation_window_lower_offset + isolation_window_upper_offset
            if full_width not in stats['isolation_window_full_widths']:
                stats['isolation_window_full_widths'][full_width] = 0
            stats['isolation_window_full_widths'][full_width] += 1

        isolation_window_lower_bound = round(isolation_window_target_mz - isolation_window_lower_offset, 2)
        isolation_window_upper_bound = round(isolation_window_target_mz + isolation_window_lower_offset, 2)
        isolation_window_str = f"{isolation_window_lower_bound}-{isolation_window_upper_bound}"

        if isolation_window_str in stats['isolation_window_last_scans']:
            scan_delta = scan_number - stats['isolation_window_last_scans'][isolation_window_str]['last_scan_number']
            if scan_delta not in stats['isolation_window_last_scans'][isolation_window_str]['scan_deltas']:
                stats['isolation_window_last_scans'][isolation_window_str]['scan_deltas'][scan_delta] = 0
            stats['isolation_window_last_scans'][isolation_window_str]['scan_deltas'][scan_delta] += 1
            stats['isolation_window_last_scans'][isolation_window_str]['last_scan_number'] = scan_number

        else:
            stats['isolation_window_last_scans'][isolation_window_str] = { 'last_scan_number': scan_number, 'scan_deltas': {} }

        try:
            precursor_ion = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']
            start_scan_time = spectrum['scanList']['scan'][0]['scan start time'] # Scan_start_time is assumed to be in minutes since pyteomics documentation says minutes
        except Exception as e:
            return

        mz_int = int(precursor_ion)

        if mz_int not in ms_one_tolerance_dict:
            ms_one_tolerance_dict[mz_int] = {"mz":[], "time":[]}

        ms_one_tolerance_dict[mz_int]['mz'].append(precursor_ion)
        ms_one_tolerance_dict[mz_int]['time'].append(start_scan_time)



    ####################################################################################################
    #### Parse the filter string
    def parse_filter_string(self,filter_string,stats):
        
        # If the filter_string is None or empty, then update a few stats and return
        if filter_string is None or filter_string == '':
            stats['high_accuracy_precursors'] = 'unknown'
            stats['fragmentation_type'] = 'unknown_fragmentation_type'
            stats['fragmentation_tag'] = 'unknown fragmentation type'
            return
        

        fragmentation_tag = '??'

        mass_accuracy = '??'
        match = re.search(r'^(\S+)',filter_string)
        if match.group(0) == 'FTMS' or match.group(0) == 'ASTMS':
            mass_accuracy = 'HR'
        elif match.group(0) == 'ITMS':
            mass_accuracy = 'LR'
        else:
            self.log_event('ERROR','UnrecognizedDetector',f"Unrecognized detector '{match.group(0)}' in filter string {filter_string}")
            return

        ms_level = 0
        match = re.search(r' ms ',filter_string)
        if match:
            ms_level = 1
        else:
            match = re.search(r' ms(\d) ',filter_string)
            if match:
                ms_level = int(match.group(1))
            else:
                self.log_event('ERROR','UnrecognizedMSLevel',f"Unrecognized MS level in filter string {filter_string}")
                return

        # See if there's fragmentation of a known type
        have_cid = 0
        have_etd = 0
        have_hcd = 0
        have_fragmentation = 0
        match = re.search(r'@hcd',filter_string)
        if match:
            have_hcd = 1
        match = re.search(r'@cid',filter_string)
        if match:
            have_cid = 1
        match = re.search(r'@etd',filter_string)
        if match:
            have_etd = 1
        match = re.search(r'@',filter_string)
        if match:
            have_fragmentation = 1
        dissociation_sum = have_cid + have_etd + have_hcd
        if have_fragmentation > dissociation_sum:
            self.log_event('ERROR','UnrecognizedFragmentation',f"Unrecognized string after @ in filter string '{filter_string}'")
            return

        #### If this is an MS0 scan, document the error and stop analyzing the spectra
        if ms_level == 0:
            if 'n_ms0_spectra' not in stats:
                stats['n_ms0_spectra'] = 0
            stats['n_ms0_spectra'] += 1
            self.log_event('ERROR','UnknownMSLevel',f"Unable to determine MS level in filter string '{filter_string}'")
            return

        #### If this is an MS1 scan, learn what we can from it
        if ms_level == 1:
            stats['n_ms1_spectra'] += 1
            spectrum_type = 'MS1'
            if mass_accuracy == 'HR':
                if stats['high_accuracy_precursors'] == 'unknown':
                    stats['high_accuracy_precursors'] = 'true'
                elif stats['high_accuracy_precursors'] == 'true' or stats['high_accuracy_precursors'] == 'multiple':
                    pass
                elif stats['high_accuracy_precursors'] == 'false':
                    stats['high_accuracy_precursors'] = 'multiple'
                else:
                    self.log_event('ERROR','InternalStateError',f"high_accuracy_precursors has a strange state '{stats['high_accuracy_precursors']}'")
                    return
            elif mass_accuracy == 'LR':
                if stats['high_accuracy_precursors'] == 'unknown':
                    stats['high_accuracy_precursors'] = 'false'
                elif stats['high_accuracy_precursors'] == 'false' or stats['high_accuracy_precursors'] == 'multiple':
                    pass
                elif stats['high_accuracy_precursors'] == 'true':
                    stats['high_accuracy_precursors'] = 'multiple'
                else:
                    self.log_event('ERROR','InternalStateError',f"high_accuracy_precursors has a strange state '{stats['high_accuracy_precursors']}'")
                    return
            else:
                self.log_event('ERROR','UnknownPrecursorScanType',f"Unknown precursor scan type '{match.group(0)}'")
                return

        #### Else it's an MSn (n>1) scan so learn what we can from that
        elif ms_level == 2:
            dissociation_sum = have_cid + have_etd + have_hcd
            stats['n_ms2_spectra'] += 1
            if dissociation_sum == 1:
                if have_hcd:
                    stats[f"n_{mass_accuracy}_HCD_spectra"] += 1
                    spectrum_type = f"{mass_accuracy}_HCD"
                    fragmentation_tag = f"{mass_accuracy} HCD"
                elif have_cid:
                    stats[f"n_{mass_accuracy}_IT_CID_spectra"] += 1
                    spectrum_type = f"{mass_accuracy}_IT_CID"
                    fragmentation_tag = f"{mass_accuracy} IT CID"
                elif have_etd:
                    stats[f"n_{mass_accuracy}_IT_ETD_spectra"] += 1
                    spectrum_type = f"{mass_accuracy}_IT_ETD"
                    fragmentation_tag = f"{mass_accuracy} IT ETD"
                else:
                    self.log_event('ERROR','UnrecognizedFragmentation',f"Did not find @hcd, @cid, @etd in filter string '{filter_string}'")
                    return

            elif dissociation_sum == 2:
                if have_etd and have_hcd:
                    stats[f"n_{mass_accuracy}_EThcD_spectra"] += 1
                    spectrum_type = f"{mass_accuracy}_EThcD"
                    fragmentation_tag = f"{mass_accuracy} EThcD"
                elif have_etd and have_cid:
                    stats[f"n_{mass_accuracy}_ETciD_spectra"] += 1
                    spectrum_type = f"{mass_accuracy}_ETciD"
                    fragmentation_tag = f"{mass_accuracy} ETciD"
                else:
                    self.log_event('ERROR','UnrecognizedFragmentation',f"Did not find known fragmentation combination in filter string '{filter_string}'")
                    return

            else:
                self.log_event('ERROR','CantParseFilterStr',f"Unable to determine dissociation type in filter string '{filter_string}'")
                return

            #### Record the fragmentation type for this run
            if spectrum_type != stats['fragmentation_type']:
                if stats['fragmentation_type'] == 'unknown':
                    stats['fragmentation_type'] = spectrum_type
                elif stats['fragmentation_type'] == 'multiple':
                    pass
                else:
                    stats['fragmentation_type'] = 'multiple'
                    self.log_event('ERROR','MultipleFragTypes',"There are multiple fragmentation types in this MS run. Split them.")


        # If we have ms_level > 2, then make a note, but this is uncharted waters
        else:
            if 'n_ms3+_spectra' not in stats:
                stats['n_ms3+_spectra'] = 0
            stats['n_ms3+_spectra'] += 1

        # Store the fragmentation_tag for use by the caller
        stats['fragmentation_tag'] = fragmentation_tag



    ####################################################################################################
    #### Add spectrum
    #def add_spectrum(self,spectrum,spectrum_type,precursor_mz):
    def add_spectrum(self,params):

        stats, precursor_mz, peaklist = params

        spectrum_type = stats['fragmentation_tag'].replace(" ", "_")

        if 'n_spectra_examined' not in self.metadata['files'][self.mzml_file]['spectra_stats']:
            self.metadata['files'][self.mzml_file]['spectra_stats']['n_spectra_examined'] = 0
        self.metadata['files'][self.mzml_file]['spectra_stats']['n_spectra_examined'] += 1

        composite_types = ["lowend", "precursor_loss"]
        
        for composite_type in composite_types:
            destination = f"{composite_type}_{spectrum_type}"

            if destination not in self.composite:
                if spectrum_type not in self.composite_spectrum_attributes[composite_type]:
                    message = f"WARNING: Not currently able to compute some metrics for spectrum type {spectrum_type}"
                    if message not in self.skipped_spectrum_types:
                        self.skipped_spectrum_types[message] = True
                        eprint(message)
                    return
                
                #print(f"INFO: Creating a composite spectrum {destination}")
                minimum = self.composite_spectrum_attributes[composite_type][spectrum_type]['minimum']
                maximum = self.composite_spectrum_attributes[composite_type][spectrum_type]['maximum']
                binsize = self.composite_spectrum_attributes[composite_type][spectrum_type]['binsize']

                array_size = int( (maximum - minimum ) / binsize ) + 1
                self.composite[destination] = { 'minimum': minimum, 'maximum': maximum, 'binsize': binsize }
                self.composite[destination]['intensities'] = numpy.zeros(array_size,dtype=numpy.float32)
                self.composite[destination]['n_peaks'] = numpy.zeros(array_size,dtype=numpy.int32)


            #### Convert the peak lists to numpy arrays
            intensity_array = numpy.asarray(peaklist['intensity array'])
            mz_array = numpy.asarray(peaklist['m/z array'])

            if composite_type == 'precursor_loss':
                mz_array = precursor_mz - mz_array


            #### Limit the numpy arrays to the region we're interested in
            intensity_array = intensity_array[ ( mz_array > self.composite[destination]['minimum'] ) & ( mz_array < self.composite[destination]['maximum'] ) ]
            mz_array = mz_array[ ( mz_array > self.composite[destination]['minimum'] ) & ( mz_array < self.composite[destination]['maximum'] ) ]


            #### Compute their bin locations and store n_peaks and intensities
            bin_array = (( mz_array - self.composite[destination]['minimum'] ) / self.composite[destination]['binsize']).astype(int)
            self.composite[destination]['n_peaks'][bin_array] += 1
            self.composite[destination]['intensities'][bin_array] += intensity_array
            

        #### The stuff below was trying to compute neutral losses. It works, but it's too slow. Need a Cython implementation I guess

        # sorted_mz_intensity = sorted(mz_intensity, key=lambda pair: pair[2], reverse=True)
        # length = 30
        # if len(mz_intensity) < 30:
            # length = len(mz_intensity)
        # clipped_mz_intensity = sorted_mz_intensity[0:(length-1)]
        # mz_intensity = sorted(clipped_mz_intensity, key=lambda pair: pair[1], reverse=False)
        # i = 0
        # while i < length-1:
            # mz_intensity[i][0] = i
            # i += 1
        # sorted_mz_intensity = sorted(mz_intensity, key=lambda pair: pair[2], reverse=True)


        # imajor_peak = 0
        # while imajor_peak < 20 and imajor_peak < len(mz_intensity):
            # ipeak,mz,major_intensity = sorted_mz_intensity[imajor_peak]
            # #print(ipeak,mz,major_intensity)
            # if major_intensity == 0.0:
                # major_intensity = 1.0
            # imz = mz
            # while ipeak > 0 and mz - imz < self.composite[destination2]['maximum'] and imz > 132 and abs(imz - precursor_mz) > 5:
                # bin = int( ( mz - imz ) / self.composite[destination2]['binsize'])
                # #print("  - ",ipeak,imz,bin,intensity)
                # self.composite[destination2]['n_peaks'][bin] += 1
                # self.composite[destination2]['intensities'][bin] += intensity
                # ipeak -= 1
                # imz = mz_intensity[ipeak][1]
                # intensity = mz_intensity[ipeak][2]
            # imajor_peak += 1


    ####################################################################################################
    #### Read header
    #### Open the mzML file and read line by line into a text buffer that we will XML parse
    def read_header(self):

        #### If the mzML is gzipped, then open with zlib, else a plain open
        match = re.search(r'\.gz$',self.mzml_file)
        if match:
            infile = gzip.open(self.mzml_file)
        else:
            infile = open(self.mzml_file)

        #### set up a text buffer to hold the mzML header
        buffer = ''
        counter = 0

        #### Keep a dict of random recognized things
        recognized_things = {}
        self.metadata['files'][self.mzml_file]['start_timestamp'] = None

        #### Read line by line
        for line in infile:

            if not isinstance(line, str):
                line = str(line, 'utf-8', 'ignore')

            #### Completely skip the <indexedmzML> tag if present
            if '<indexedmzML ' in line:
                continue

            #### Look for first tag after the header and end when found
            if '<run ' in line:
                #### Extract the startTimeStamp
                match = re.search(r'startTimeStamp="(.+?)"', line)
                if match:
                    self.metadata['files'][self.mzml_file]['start_timestamp'] = match.group(1)
                break

            #### Look for some known artificts that are useful to record
            if 'softwareRef="tdf2mzml"' in line:
                recognized_things['converted with tdf2mzml'] = True
                recognized_things['likely a timsTOF'] = True

            #### Update the buffer and counter
            if counter > 0:
                buffer += line
            counter += 1

        #### Close file and process the XML in the buffer
        infile.close()

        #### Finish the XML by closing the tags
        buffer += '  </mzML>\n'

        #### Get the root of the XML and the namespace
        xmlroot = etree.fromstring(buffer)
        namespace = xmlroot.nsmap
        #print(namespace)
        if None in namespace:
            namespace = '{'+namespace[None]+'}'
        else:
            namespace = ''

        #### Create a reference of instruments we know
        instrument_by_category = {
            'pureHCD': [
                'MS:1000649|Exactive',
                'MS:1002526|Exactive Plus',
                'MS:1001911|Q Exactive',
                'MS:1002993|Q Exactive Focus',
                'MS:1002523|Q Exactive HF',
                'MS:1002877|Q Exactive HF-X',
                'MS:1002634|Q Exactive Plus',
                'MS:1003245|Q Exactive UHMR',
                'MS:1003378|Orbitrap Astral',
                'MS:1003442|Orbitrap Astral Zoom'
            ],
            'ion_trap': [
                'MS:1000447|LTQ',
                'MS:1000855|LTQ Velos',
                'MS:1000856|LTQ Velos/ETD',
                'MS:1000854|LTQ XL',
                'MS:1000638|LTQ XL ETD',
                'MS:1000167|LCQ Advantage',
                'MS:1000168|LCQ Classic',
                'MS:1000169|LCQ Deca XP Plus',
                'MS:1000554|LCQ Deca',
                'MS:1000578|LCQ Fleet',
                'MS:1001909|Velos Plus',
                'MS:1003495|Velos Pro'
            ],
            'variable': [ 
                'MS:1000448|LTQ FT', 
                'MS:1000557|LTQ FT Ultra', 
                'MS:1000449|LTQ Orbitrap',
                'MS:1002835|LTQ Orbitrap Classic', 
                'MS:1000555|LTQ Orbitrap Discovery', 
                'MS:1001742|LTQ Orbitrap Velos', 
                'MS:1003499|LTQ Orbitrap Velos/ETD', 
                'MS:1000556|LTQ Orbitrap XL',
                'MS:1000639|LTQ Orbitrap XL ETD',
                'MS:1003356|Orbitrap Ascend',
                'MS:1003029|Orbitrap Eclipse',
                'MS:1001910|Orbitrap Elite', 
                'MS:1002994|Orbitrap Excedion Pro', 
                'MS:1003095|Orbitrap Exploris 120',
                'MS:1003094|Orbitrap Exploris 240',
                'MS:1003028|Orbitrap Exploris 480',
                'MS:1003423|Orbitrap Exploris GC 240',
                'MS:1002992|Orbitrap Exploris GC-MS',
                'MS:1002416|Orbitrap Fusion', 
                'MS:1002417|Orbitrap Fusion ETD', 
                'MS:1002732|Orbitrap Fusion Lumos',
                'MS:1003096|Orbitrap Velos Pro', 
                'MS:1000031|instrument model',
                'MS:1000483|Thermo Fisher Scientific instrument model'
            ],
            'QTOF': [
                'MS:1000126|Waters instrument model',
                'MS:1000122|Bruker Daltonics instrument model',
                'MS:1001547|Bruker Daltonics maXis series',
                'MS:1003123|Bruker Daltonics timsTOF series',
                'MS:1000121|AB SCIEX instrument model',
                'MS:1002583|TripleTOF 4600',
                'MS:1000932|TripleTOF 5600',
                'MS:1002584|TripleTOF 5600+',
                'MS:1002533|TripleTOF 6600',
                'MS:1003456|TripleTOF 6600+',
                'MS:1003293|ZenoTOF 7600',
                'MS:1003444|ZenoTOF 7600+',
                'MS:1003443|ZenoTOF 8600'
            ]
        }

        #### Some special cases to look out for
        special_cases = {
            'MS:1000566': {
                'exact_name': 'ISB mzXML format',
                'rank': '5',
                'category': 'unknown',
                'final_name': 'Conversion to mzML from mzXML from unknown instrument' },
            'MS:1000562': {
                'exact_name': 'ABI WIFF format',
                'rank': '2',
                'category': 'QTOF',
                'final_name': 'SCIEX instrument model' },
            'MS:1000551': {
                'exact_name': 'Analyst',
                'rank': '3',
                'category': 'QTOF',
                'final_name': 'SCIEX instrument model' }
        }

        #### Restructure it into a dict by PSI-MS identifier
        instrument_attributes = {}
        for instrument_category in instrument_by_category:
            for instrument_string in instrument_by_category[instrument_category]:
                accession,name = instrument_string.split('|')
                instrument_attributes[accession] = { 'category': instrument_category, 'accession': accession, 'name': name }

        #### Get all the CV params in the header and look for ones we know about
        cv_params = xmlroot.findall(f'.//{namespace}cvParam')
        found_instrument = 0
        model_data = None
        for cv_param in cv_params:
            accession = cv_param.get('accession')
            if accession in instrument_attributes:
                #### Store the attributes about this instrument model
                model_data = {
                    'accession': accession,
                    'name': instrument_attributes[accession]['name'],
                    'category': instrument_attributes[accession]['category'],
                    'rank': '1'
                }
                if 'likely a timsTOF' in recognized_things:
                    model_data['inferred_name'] = 'timsTOF'
                self.metadata['files'][self.mzml_file]['instrument_model'] = model_data
                found_instrument = 1

            #### Handle some special cases
            if accession in special_cases:
                if model_data is None or special_cases[accession]['rank'] < model_data['rank']:
                    model_data = {
                        'accession': accession,
                        'name': special_cases[accession]['final_name'],
                        'category': special_cases[accession]['category'],
                        'rank': special_cases[accession]['rank']
                    }
                    self.metadata['files'][self.mzml_file]['instrument_model'] = model_data
                    found_instrument = 1

        #### If none are an instrument we know about about, ask for help
        if not found_instrument:
            self.log_event('WARNING','UnrecogInstr',"Did not find the instrument CV term. Please update RunAssessor with information about this instrument. Proceding as if it were a QTOF type instrument.")
            model_data = {
                'accession': None,
                'name': 'unknown',
                'category': 'QTOF'
            }
            self.metadata['files'][self.mzml_file]['instrument_model'] = model_data
            return


        #### Read and store the referenceable param groups
        referenceable_param_group_list = xmlroot.findall(f'.//{namespace}referenceableParamGroup')
        for referenceable_param_group in referenceable_param_group_list:
            group_id = referenceable_param_group.get('id')
            self.referenceable_param_group_list[group_id] = {}
            for subelement in list(referenceable_param_group):
                name = subelement.get('name')
                value = subelement.get('value')
                self.referenceable_param_group_list[group_id][name] = value

        return model_data



    ####################################################################################################
    #### Assess low end composite spectra
    def assess_lowend_composite_spectra(self):
        #with open('zzcomposite_spectrum.pickle','rb') as infile:
        #    spec = pickle.load(infile)
        spec = self.composite

        ### File path for table including additional peaks
        additionalPeaks = os.path.dirname(os.path.abspath(__file__)) + "/Additional_Peaks.tsv"

        ###Dictionary to hold upper and lower sigma_ppm relative to theoretical value
        self.all_3sigma_values_away = {"Status": False,"upper_ppm": [], "lower_ppm": [], "upper_mz":[], "lower_mz":[]}

        supported_composite_type_list = [f"lowend_{composite}" for composite in self.supported_composite_type_list]

        composite_type_list = []
        for supported_composite_type in supported_composite_type_list:
            if supported_composite_type in spec:
                
                composite_type_list.append(supported_composite_type)

        if not composite_type_list:
            self.log_event('WARNING','UnknownFragmentation',"Not able to determine more about these spectra yet.")
            return
        
        for composite_type in composite_type_list:
            minimum = spec[composite_type]['minimum']
            maximum = spec[composite_type]['maximum']
            binsize = spec[composite_type]['binsize']
            spec[composite_type]['mz'] = numpy.arange(minimum, maximum + binsize, binsize)
            self.metadata['files'][self.mzml_file].setdefault('lowend_peaks',{})

            ROIs = {
                'TMT6_126': { 'type': 'TMT', 'mz': 126.127725, 'initial_window': 0.01 },
                'TMT6_127': { 'type': 'TMT', 'mz': 127.124760, 'initial_window': 0.01 },
                'TMT6_128': { 'type': 'TMT', 'mz': 128.134433, 'initial_window': 0.01 },
                'TMT6_129': { 'type': 'TMT', 'mz': 129.131468, 'initial_window': 0.01 },
                'TMT6_130': { 'type': 'TMT', 'mz': 130.141141, 'initial_window': 0.01 },
                'TMT6_131': { 'type': 'TMT', 'mz': 131.138176, 'initial_window': 0.01 },
                'TMT10_127_C': {'type': 'TMT10', 'mz': 127.131081, 'initial_window': 0.01 },
                'TMT10_128_N': {'type': 'TMT10', 'mz': 128.128116, 'initial_window': 0.01 },
                'TMT10_129_C': {'type': 'TMT10', 'mz': 129.137790, 'initial_window': 0.01 },
                'TMT10_130_N': {'type': 'TMT10', 'mz': 130.134825, 'initial_window': 0.01 },
                'TMT11_131_C': {'type': 'TMT10', 'mz': 131.144499, 'initial_window': 0.01 },
                'TMT6plex': { 'type': 'TMT6', 'mz': 230.1702, 'initial_window': 0.01 },
                'TMT6plex+H2O': { 'type': 'TMT6', 'mz': 248.18079, 'initial_window': 0.01 },
                'TMTpro': { 'type': 'TMTpro', 'mz': 305.2144, 'initial_window': 0.01 },
                'TMTpro+H2O': { 'type': 'TMTpro', 'mz': 323.22499, 'initial_window': 0.01 },
                'iTRAQ_114': { 'type': 'iTRAQ', 'mz': 114.11068, 'initial_window': 0.01 },
                'iTRAQ_115': { 'type': 'iTRAQ', 'mz': 115.107715, 'initial_window': 0.01 },
                'iTRAQ_116': { 'type': 'iTRAQ', 'mz': 116.111069, 'initial_window': 0.01 },
                'iTRAQ_117': { 'type': 'iTRAQ', 'mz': 117.114424, 'initial_window': 0.01 },
                'iTRAQ4plex': { 'type': 'iTRAQ4', 'mz': 145.109, 'initial_window': 0.01 },
                'iTRAQ8_118': { 'type': 'iTRAQ8', 'mz': 118.111459, 'initial_window': 0.01 },
                'iTRAQ8_119': { 'type': 'iTRAQ8', 'mz': 119.114814, 'initial_window': 0.01 },
                'iTRAQ8_121': { 'type': 'iTRAQ8', 'mz': 121.121524, 'initial_window': 0.01 },
                'iTRAQ8_113': { 'type': 'iTRAQ8', 'mz': 113.107325, 'initial_window': 0.01 },
                'iTRAQ8plex': { 'type': 'iTRAQ8', 'mz': 304.205360+1.00727, 'initial_window': 0.01 }
            }
            
            #Read additional peaks into ROI dictionary 
            # 0: mean monoisotopic m/z       
            # 1: average intensity           
            # 2: # MS runs                  
            # 3: percentage of spectra       
            # 4: theoretical m/z             
            # 5: delta PPM                  
            # 6: primary identification      
            # 7: other identifications    
            with open(additionalPeaks, "r") as f:
                next(f)
                for line in f:
                    parts = line.strip().split("\t")
                    if r"a{AI}" in parts[6] or r"b{AK}-NH3" in parts[6]: # On 7/11/2025 we decided to not take these ions (m/z might be wrong)
                        pass
                    elif float(parts[4]) > 100:
                        ROIs[parts[6]] = {'type':'fragment_ion', "mz":float(parts[4]), 'initial_window': 0.01}
            
            #### What should we look at for ion trap data?
            if composite_type in [ 'lowend_LR_IT_CID', 'lowend_LR_HCD' ]:
                ROIs = {
                    'TMT6plex': { 'type': 'TMT6', 'mz': 230.1702, 'initial_window': 1.0 },
                    'TMT6plex+H2O': { 'type': 'TMT6', 'mz': 248.18079, 'initial_window': 1.0 },
                    'TMT6_y1K': { 'type': 'TMT6', 'mz': 376.2757, 'initial_window': 1.0 },
                    'TMTpro': { 'type': 'TMTpro', 'mz': 305.2144, 'initial_window': 1.0 },
                    'TMTpro+H2O': { 'type': 'TMTpro', 'mz': 323.22499, 'initial_window': 1.0 },
                    'iTRAQ4_y1K': { 'type': 'iTRAQ4', 'mz': 376.2757, 'initial_window': 1.0 },
                    'iTRAQ8_y1K': { 'type': 'iTRAQ8', 'mz': 451.3118, 'initial_window': 1.0 }
                }
                with open(additionalPeaks, "r") as f:
                    next(f)
                    for line in f:
                        parts = line.strip().split("\t")
                        if r"a{AI}" in parts[6] or r"b{AK}-NH3" in parts[6]: # On 7/11/2025 we decided to not take these ions (m/z might be wrong)
                            pass
                        elif float(parts[4]) > 200:
                            ROIs[parts[6]] = {'type':'fragment_ion', "mz":float(parts[4]), 'initial_window': 1.0}


            for ROI in ROIs:

                #### Set thresholds for peak finding
                if ROI.startswith("TMT") or ROI.startswith("iTRAQ"):
                    thresholds = {
                        'min_spectra': 100,
                        'min_percentage': 0.1
                    }
                else:
                    thresholds = {
                        'min_spectra': 100,
                        'min_percentage': 0.01
                    }

                #### Look for the largest peak
                peak = self.find_peak(spec,ROIs,ROI,composite_type,thresholds)

                ROIs[ROI]['peak'] = peak

                try:
                    ### Finds the number of standard deviations away the peak is from the theoretical mz value
                    sigma_away_from_theoretical = (ROIs[ROI]['peak']['fit']["delta_ppm"])/(ROIs[ROI]['peak']['fit']["sigma_ppm"])
                    ### Assigning the information in the dictionary
                    ROIs[ROI]['peak']['assessment']['# of sigma_ppm(s) away from theoretical peak (ppm)'] = sigma_away_from_theoretical
                    
                except:
                    pass

                #If a 3 standard deviation is availible from the guassian curve, check to see if its distance is a max or min in the whole spectra set

                try:
                    if self.fragmentation_units in {"mz", "ppm"}:
                        units = self.fragmentation_units
                    else:
                        if composite_type == "lowend_HR_HCD":
                            units = "ppm"
                        elif composite_type in {"lowend_LR_IT_CID"}:
                            units = "mz"
                        else:
                            units = None

                    if units:
                        upper = peak['extended'][f'three_sigma_{units}_upper']
                        lower = peak['extended'][f'three_sigma_{units}_lower']

                        self.all_3sigma_values_away[f'upper_{units}'].append(upper)
                        self.all_3sigma_values_away[f'lower_{units}'].append(lower)
                        self.all_3sigma_values_away['Status'] = True
                                    

                except KeyError:
                    pass

            
            self.metadata['files'][self.mzml_file]['lowend_peaks'][composite_type] = ROIs
        return ROIs
    

    ####################################################################################################
    #### Assess neutral loss composite spectra
    def assess_neutral_loss_composite_spectra(self):
        #with open('zzcomposite_spectrum.pickle','rb') as infile:
        #    spec = pickle.load(infile)
        spec = self.composite

        supported_composite_type_list = [f"precursor_loss_{composite}" for composite in self.supported_composite_type_list]

        self.metadata['files'][self.mzml_file].setdefault('neutral_loss_peaks', {})
        composite_type_list = []
        for supported_composite_type in supported_composite_type_list:
            if supported_composite_type in spec:
                
                composite_type_list.append(supported_composite_type)

        if not composite_type_list:
            self.log_event('WARNING','UnknownFragmentation',"Not able to determine more about these spectra yet.")
            return
        for composite_type in composite_type_list:
            minimum = spec[composite_type]['minimum']
            maximum = spec[composite_type]['maximum']
            binsize = spec[composite_type]['binsize']
            spec[composite_type]['mz'] = numpy.arange(minimum, maximum + binsize, binsize)

            singly_charged_peaks = [
                { 'name': 'phospho', 'type': 'phospho', 'mass': 79.966331 },
                { 'name': 'water', 'type': 'water_loss', 'mass': 18.0105647 },
                { 'name': 'phosphoric_acid', 'type': 'phosphoric_acid_loss', 'mass': 97.9768957 }
                #{ 'name': 'ammonia', 'type': 'ammonia_loss', 'mass': 17.026549105 }
            ]
            ROIs = {}
            charges = [1, 2, 3]
            for item in singly_charged_peaks:
                for z in charges:
                    ROI_name = f"{item['name']}_z{z}"
                    mz = (item['mass'] / z)
                    type = f"{item['type']}"

                    ROIs[ROI_name] = {'type': type, 'mz': mz, 'initial_window': 0.01}

            #### What should we look at for ion trap data?
            if composite_type in [ 'precursor_loss_LR_IT_CID', 'precursor_loss_LR_HCD' ]:
                ROIs = {}
                charges = [1, 2, 3]
                for item in singly_charged_peaks:
                    for z in charges:
                        ROI_name = f"{item['name']}_z{z}"
                        mz = (item['mass'] / z)
                        type = f"{item['type']}"

                        ROIs[ROI_name] = {'type': type, 'mz': mz, 'initial_window': 1.0}

            #### Set thresholds for peak finding
            thresholds = {
                'min_spectra': 100,
                'min_percentage': 0.1
            }
            if composite_type == 'precursor_loss_HR_HCD':
                thresholds['min_percentage'] = 0


            for ROI in ROIs:

                #### Look for the largest peak
                peak = self.find_peak(spec,ROIs,ROI,composite_type,thresholds)
                ROIs[ROI]['peak'] = peak

                try:
                    del(peak['fit']['sigma_ppm'])
                    del(peak['fit']['delta_ppm'])
                except:
                    pass

            
            self.metadata['files'][self.mzml_file]['neutral_loss_peaks'][composite_type] = ROIs
        return ROIs


    ####################################################################################################
    #### Peak finding routine
    def find_peak(self,spec,ROIs,ROI,composite_type,thresholds):
        minimum = spec[composite_type]['minimum']
        #maximum = spec[composite_type]['maximum']
        binsize = spec[composite_type]['binsize']
        center = ROIs[ROI]['mz']
        range = ROIs[ROI]['initial_window']
        first_bin = int((center - range/2.0 - minimum) / binsize)
        last_bin = int((center + range/2.0 - minimum) / binsize)

        peak = { 'mode_bin': { 'ibin': 0, 'mz': 0, 'intensity': 0, 'n_spectra': 0 } }

        #### Find the tallest peak in the window
        ibin = first_bin
        while ibin < last_bin:
            if spec[composite_type]['intensities'][ibin] > peak['mode_bin']['intensity']:
                peak['mode_bin']['ibin'] = ibin
                peak['mode_bin']['mz'] = float(spec[composite_type]['mz'][ibin])
                peak['mode_bin']['intensity'] = float(spec[composite_type]['intensities'][ibin])
                peak['mode_bin']['n_spectra'] = int(spec[composite_type]['n_peaks'][ibin])
            ibin += 1

        #### If there were 0 peaks found, return a verdict now
        if peak['mode_bin']['ibin'] == 0:
            peak['assessment'] = { 'is_found': False, 'fraction': 0.0 }
            return(peak)

        #### Grow the extent to get most of the signal
        ibin = peak['mode_bin']['ibin']
        iextent = 1
        done = 0
        prev_n_spectra = peak['mode_bin']['n_spectra']
        peak['extended'] = { 'extent': 0, 'intensity': float(peak['mode_bin']['intensity']), 'n_spectra': int(peak['mode_bin']['n_spectra']) }
        #print(spec[composite_type]['n_peaks'][ibin-10:ibin+10])
        while not done:
            prev_intensity = peak['extended']['intensity']
            prev_n_spectra = peak['extended']['n_spectra']
            new_intensity = prev_intensity + spec[composite_type]['intensities'][ibin-iextent] + spec[composite_type]['intensities'][ibin+iextent]
            new_n_spectra = prev_n_spectra + spec[composite_type]['n_peaks'][ibin-iextent] + spec[composite_type]['n_peaks'][ibin+iextent]
            #print(new_intensity/prev_intensity, new_n_spectra/prev_n_spectra)
            peak['extended'] = { 'extent': iextent, 'intensity': float(new_intensity), 'n_spectra': int(new_n_spectra) }
            if new_n_spectra/prev_n_spectra < 1.02:
                #print('Extent reached by n_spectra')
                done = 1
            if iextent > 15:
                #print('Extent reached by max extent')
                done = 1
            iextent += 1


        min_spectra = thresholds['min_spectra']
        min_percentage = thresholds['min_percentage']

        if peak['extended']['n_spectra'] > min_spectra and peak['extended']['n_spectra'] >= self.metadata['files'][self.mzml_file]['spectra_stats']['n_spectra'] * min_percentage:
            extent = peak['extended']['extent'] * 2
            x = spec[composite_type]['mz'][ibin-extent:ibin+extent]
            y = spec[composite_type]['intensities'][ibin-extent:ibin+extent]

            n = len(x)
            center = int(n/2)
            binsize = x[center]-x[center-1]

            try:
            #if 1:
                popt,pcov = curve_fit(gaussian_function,x,y,p0=[y[center],x[center],binsize, 0.0])
                sigma_ppm = popt[2]/ROIs[ROI]['mz']*1e6
                
                peak['fit'] = { 'mz': popt[1], 'intensity': popt[0], 'sigma_mz': popt[2], 'delta_mz': popt[1]-ROIs[ROI]['mz'], 'delta_ppm': (popt[1]-ROIs[ROI]['mz'])/ROIs[ROI]['mz']*1e6, 'sigma_ppm': sigma_ppm, "y_offset": popt[3]}
                if composite_type == "lowend_HR_HCD":
                    peak['extended']['three_sigma_ppm_lower'] = peak['fit']['delta_ppm']-sigma_ppm*3
                    peak['extended']['three_sigma_ppm_upper'] = peak['fit']['delta_ppm']+3*sigma_ppm
                    peak['extended']['three_sigma_mz_lower'] = peak['fit']['delta_mz']-popt[2]*3
                    peak['extended']['three_sigma_mz_upper'] = peak['fit']['delta_mz']+popt[2]*3

                elif composite_type in [ 'lowend_LR_IT_CID', 'lowend_LR_HCD' ]:
                    peak['extended']['three_sigma_mz_lower'] = peak['fit']['delta_mz']-popt[2]*3
                    peak['extended']['three_sigma_mz_upper'] = peak['fit']['delta_mz']+popt[2]*3
                    peak['extended']['three_sigma_ppm_lower'] = peak['fit']['delta_ppm']-sigma_ppm*3
                    peak['extended']['three_sigma_ppm_upper'] = peak['fit']['delta_ppm']+3*sigma_ppm
                
                elif composite_type.startswith("precursor_"):
                    pass
                
                peak['assessment'] = { 'is_found': True, 'fraction': 0.0, 'comment': 'Peak found and fit' }
            except:
            #else:
                peak['assessment'] = { 'is_found': False, 'fraction': 0.0, 'comment': 'Gaussian fit failed to converge' }

            if 0:
                plt.step(x+binsize/2.0,y,'b+:')
                plt.plot(x,gaussian_function(x,*popt),'ro:')
                title = "Fit results: mz = %.4f, sigma = %.4f" % (popt[1], popt[2])
                plt.title(title)
                plt.show()

        else:
            peak['assessment'] = { 'is_found': False, 'fraction': 0.0, 'comment': 'Too few peaks found in ROI' }

        

        return peak


    ####################################################################################################
    #### Assess Regions of Interest to make calls
    def assess_ROIs(self):

        possible_frag_type = ['n_HR_HCD_spectra', 'n_LR_HCD_spectra', 'n_HR_IT_CID_spectra', 'n_LR_IT_CID_spectra', 'n_HR_IT_ETD_spectra', 'n_LR_IT_ETD_spectra', 'n_HR_EThcD_spectra', 'n_LR_EThcD_spectra', 'n_HR_ETciD_spectra', 'n_LR_ETciD_spectra', 'n_HR_QTOF_spectra']
        fragmentation_types = []
        for frag_type in self.metadata['files'][self.mzml_file]['spectra_stats']:
            if frag_type in possible_frag_type:
                if int(self.metadata['files'][self.mzml_file]['spectra_stats'][frag_type]) > 0:
                    fragmentation_types.append(frag_type.replace("n_","").replace("_spectra", ""))
                    
  

        for fragmentations in fragmentation_types:

            results = { 'labeling': { 'scores': { 'TMT': 0, 'TMT6': 0, 'TMT10': 0, 'TMTpro': 0, 'iTRAQ': 0, 'iTRAQ4': 0, 'iTRAQ8': 0} } }
            TMT_all = 0
            iTRAQ_all = 0
            #### Determine what the denominator of MS2 spectra should be
            n_ms2_spectra = self.metadata['files'][self.mzml_file]['spectra_stats']["n_" + fragmentations + "_spectra"]
            if ( 'n_HCD_spectra' in self.metadata['files'][self.mzml_file]['spectra_stats'] and
                    self.metadata['files'][self.mzml_file]['spectra_stats']['n_HCD_spectra'] > 0 and 
                    self.metadata['files'][self.mzml_file]['spectra_stats']['n_HCD_spectra'] <
                    self.metadata['files'][self.mzml_file]['spectra_stats']['n_ms2_spectra'] ):
                n_ms2_spectra = self.metadata['files'][self.mzml_file]['spectra_stats']['n_HCD_spectra']
            
            #### Determine the fragmentation_type
            results['fragmentation_type'] = fragmentations
            
            self.metadata['files'][self.mzml_file].setdefault('summary', {})
            self.metadata['files'][self.mzml_file]['summary'][fragmentations] = results

            #### If standard deviations have been found, find the difference between them, if none have been found, mention it
            self.metadata['files'][self.mzml_file]['summary'][fragmentations].setdefault('tolerance', {})

            if self.all_3sigma_values_away['Status']:
                percentile = 90 #Picks this percentile for the three_sigma value
                try:
                    if self.fragmentation_units in {"mz", "ppm"}:
                        units = self.fragmentation_units
                    else:
                        if fragmentations.startswith('HR'):
                            units = "ppm"
                        elif fragmentations.startswith('LR'):
                            units = "mz"
                        else:
                            units = None

                    if units:
                        lower = self.all_3sigma_values_away[f'lower_{units}']
                        upper = self.all_3sigma_values_away[f'upper_{units}']
                        
                        three_sigma_upper = numpy.percentile(upper, percentile)
                        three_sigma_lower = numpy.percentile(lower, 100-percentile)

                        self.metadata['files'][self.mzml_file]['summary'][fragmentations]['tolerance'][f'fragment_tolerance_{units}_lower'] = three_sigma_lower
                        self.metadata['files'][self.mzml_file]['summary'][fragmentations]['tolerance'][f'fragment_tolerance_{units}_upper'] = three_sigma_upper
                        self.metadata['files'][self.mzml_file]['summary'][fragmentations]['tolerance']['number of peaks with fragment_tolerance'] = len(upper)
                        if units == 'ppm':
                            recommended_precursor = math.ceil(math.sqrt(self.ppm_error**2 + (max(abs(three_sigma_lower), abs(three_sigma_upper))**2)))
                        elif units == 'mz':
                            recommended_precursor = round(math.sqrt(self.mz_error**2 + (max(abs(three_sigma_lower), abs(three_sigma_upper)))**2), 4)
                            units = 'm/z'
                        self.metadata['files'][self.mzml_file]['summary'][fragmentations]['tolerance']['recommended fragment tolerance'] = recommended_precursor
                        self.metadata['files'][self.mzml_file]['summary'][fragmentations]['tolerance']['recommended fragment tolerance units'] = units

                        if len(upper) < 5:
                            self.metadata['files'][self.mzml_file]['summary'][fragmentations]['tolerance']['warning'] = "too few peaks found, overall fragment tolerance cannot be accurately computed"
        
                        elif abs(three_sigma_lower) > 20 or abs(three_sigma_upper) > 20:
                            self.metadata['files'][self.mzml_file]['summary'][fragmentations]['tolerance']['warning'] = "overall fragment tolerance levels are too high, peak fitting may not be accurate"
                        
                    

                except:
                    self.metadata['files'][self.mzml_file]['summary'][fragmentations]['tolerance']['warning'] = 'no tolerance values recorded'
            else:
                self.metadata['files'][self.mzml_file]['summary'][fragmentations]['tolerance']['warning'] = 'no tolerance values recorded'

           



            #### If no ROIs were computed, cannot continue
            if 'lowend_peaks' not in self.metadata['files'][self.mzml_file] or 'neutral_loss_peaks' not in self.metadata['files'][self.mzml_file]:
                full_results = {"fragmentation type": fragmentations, "call": "unavailable", "fragmentation tolerance": 'unavailable', "has water_loss":'unavailable', "has phospho_spectra": 'unavailable'}
                self.metadata['files'][self.mzml_file]['summary']['combined summary'] = full_results
                return

  
            ROI_type = "lowend_" + fragmentations
            for ROI in self.metadata['files'][self.mzml_file]['lowend_peaks'][ROI_type]:
                ROIobj = self.metadata['files'][self.mzml_file]['lowend_peaks'][ROI_type][ROI]
                if ROIobj['peak']['assessment']['is_found']:
                    type = ROIobj['type']
                    n_spectra = ROIobj['peak']['extended']['n_spectra']
                    ROIobj['peak']['assessment']['fraction'] = n_spectra / n_ms2_spectra
                    try:
                        results['labeling']['scores'][type] += ROIobj['peak']['assessment']['fraction']
                    except:
                        pass
            
            
            #### Detect water loss z=2 or not
            loss_type = 'precursor_loss_' + fragmentations
            if self.metadata['files'][self.mzml_file]['neutral_loss_peaks'][loss_type]['water_z2']['peak']['mode_bin']['n_spectra'] >= 50:
                if 'fit' in self.metadata['files'][self.mzml_file]['neutral_loss_peaks'][loss_type]['water_z2']['peak']:
                    if self.metadata['files'][self.mzml_file]['neutral_loss_peaks'][loss_type]['water_z2']['peak']['fit']['sigma_mz'] < 0.1:
                        self.metadata['files'][self.mzml_file]['summary'][fragmentations]['has water_loss'] = True
                    else:
                        self.metadata['files'][self.mzml_file]['summary'][fragmentations]['has water_loss'] = False
                    if ('LR' in fragmentations or 'HR_IT_ETD' in fragmentations) and abs(self.metadata['files'][self.mzml_file]['neutral_loss_peaks'][loss_type]['water_z2']['peak']['fit']['delta_mz']) > 0.4:
                        self.metadata['files'][self.mzml_file]['summary'][fragmentations]['has water_loss'] = False
                    elif 'HR' in fragmentations and abs(self.metadata['files'][self.mzml_file]['neutral_loss_peaks'][loss_type]['water_z2']['peak']['fit']['delta_mz']) > 0.008:
                        self.metadata['files'][self.mzml_file]['summary'][fragmentations]['has water_loss'] = False
                else:
                    self.metadata['files'][self.mzml_file]['summary'][fragmentations]['has water_loss'] = False
            else:
                self.metadata['files'][self.mzml_file]['summary'][fragmentations]['has water_loss'] = False


            epsilon = 1e-10
            if self.metadata['files'][self.mzml_file]['neutral_loss_peaks'][loss_type]['water_z2']['peak']['assessment']['is_found'] is True and self.metadata['files'][self.mzml_file]['neutral_loss_peaks'][loss_type]['phosphoric_acid_z2']['peak']['assessment']['is_found'] is True:
                
                try:
                    self.metadata['files'][self.mzml_file]['summary'][fragmentations]['intensity of z=2 phosphoric_acid'] = self.metadata['files'][self.mzml_file]['neutral_loss_peaks'][loss_type]['phosphoric_acid_z2']['peak']['fit']['intensity'] - self.metadata['files'][self.mzml_file]['neutral_loss_peaks'][loss_type]['phosphoric_acid_z2']['peak']['fit']['y_offset']
                    self.metadata['files'][self.mzml_file]['summary'][fragmentations]['intensity of z=2 water_loss'] = self.metadata['files'][self.mzml_file]['neutral_loss_peaks'][loss_type]['water_z2']['peak']['fit']['intensity'] - self.metadata['files'][self.mzml_file]['neutral_loss_peaks'][loss_type]['water_z2']['peak']['fit']['y_offset']
                    self.metadata['files'][self.mzml_file]['summary'][fragmentations]['z=2 phosphoric_acid to z=2 water_loss intensity ratio'] = self.metadata['files'][self.mzml_file]['summary'][fragmentations]['intensity of z=2 phosphoric_acid'] / self.metadata['files'][self.mzml_file]['summary'][fragmentations]['intensity of z=2 water_loss']
                    self.metadata['files'][self.mzml_file]['summary'][fragmentations]['absolute difference in delta m/z for z=2 phosphoric_acid_loss and z=2 water_loss'] = abs( self.metadata['files'][self.mzml_file]['neutral_loss_peaks'][loss_type]['phosphoric_acid_z2']['peak']['fit']['delta_mz'] - self.metadata['files'][self.mzml_file]['neutral_loss_peaks'][loss_type]['water_z2']['peak']['fit']['delta_mz'] )
                except:
                    self.metadata['files'][self.mzml_file]['summary'][fragmentations]['intensity of z=2 phosphoric_acid'] = 'n/a'
                    self.metadata['files'][self.mzml_file]['summary'][fragmentations]['intensity of z=2 water_loss'] = 'n/a'
                    self.metadata['files'][self.mzml_file]['summary'][fragmentations]['z=2 phosphoric_acid to z=2 water_loss intensity ratio'] = 'n/a'
                    self.metadata['files'][self.mzml_file]['summary'][fragmentations]['absolute difference in delta m/z for z=2 phosphoric_acid_loss and z=2 water_loss'] = 'n/a'
            
            elif self.metadata['files'][self.mzml_file]['neutral_loss_peaks'][loss_type]['water_z2']['peak']['assessment']['is_found'] is True and self.metadata['files'][self.mzml_file]['neutral_loss_peaks'][loss_type]['phosphoric_acid_z2']['peak']['assessment']['is_found'] is False:
                
                try:
                    self.metadata['files'][self.mzml_file]['summary'][fragmentations]['intensity of z=2 phosphoric_acid'] = epsilon
                    self.metadata['files'][self.mzml_file]['summary'][fragmentations]['intensity of z=2 water_loss'] = self.metadata['files'][self.mzml_file]['neutral_loss_peaks'][loss_type]['water_z2']['peak']['fit']['intensity'] - self.metadata['files'][self.mzml_file]['neutral_loss_peaks'][loss_type]['water_z2']['peak']['fit']['y_offset']
                    self.metadata['files'][self.mzml_file]['summary'][fragmentations]['z=2 phosphoric_acid to z=2 water_loss intensity ratio'] = self.metadata['files'][self.mzml_file]['summary'][fragmentations]['intensity of z=2 phosphoric_acid'] / self.metadata['files'][self.mzml_file]['summary'][fragmentations]['intensity of z=2 water_loss']
                    self.metadata['files'][self.mzml_file]['summary'][fragmentations]['absolute difference in delta m/z for z=2 phosphoric_acid_loss and z=2 water_loss'] = 'n/a'
                except:
                    self.metadata['files'][self.mzml_file]['summary'][fragmentations]['intensity of z=2 phosphoric_acid'] = 'n/a'
                    self.metadata['files'][self.mzml_file]['summary'][fragmentations]['intensity of z=2 water_loss'] = 'n/a'
                    self.metadata['files'][self.mzml_file]['summary'][fragmentations]['z=2 phosphoric_acid to z=2 water_loss intensity ratio'] = 'n/a'
                    self.metadata['files'][self.mzml_file]['summary'][fragmentations]['absolute difference in delta m/z for z=2 phosphoric_acid_loss and z=2 water_loss'] = 'n/a'
            
            elif self.metadata['files'][self.mzml_file]['neutral_loss_peaks'][loss_type]['water_z2']['peak']['assessment']['is_found'] is False and self.metadata['files'][self.mzml_file]['neutral_loss_peaks'][loss_type]['phosphoric_acid_z2']['peak']['assessment']['is_found'] is True:
                
                try:
                    self.metadata['files'][self.mzml_file]['summary'][fragmentations]['intensity of z=2 phosphoric_acid'] = self.metadata['files'][self.mzml_file]['neutral_loss_peaks'][loss_type]['phosphoric_acid_z2']['peak']['fit']['intensity'] - self.metadata['files'][self.mzml_file]['neutral_loss_peaks'][loss_type]['phosphoric_acid_z2']['peak']['fit']['y_offset']
                    self.metadata['files'][self.mzml_file]['summary'][fragmentations]['intensity of z=2 water_loss'] = epsilon
                    self.metadata['files'][self.mzml_file]['summary'][fragmentations]['z=2 phosphoric_acid to z=2 water_loss intensity ratio'] = self.metadata['files'][self.mzml_file]['summary'][fragmentations]['intensity of z=2 phosphoric_acid'] / self.metadata['files'][self.mzml_file]['summary'][fragmentations]['intensity of z=2 water_loss']
                    self.metadata['files'][self.mzml_file]['summary'][fragmentations]['absolute difference in delta m/z for z=2 phosphoric_acid_loss and z=2 water_loss'] = 'n/a'
                except:
                    self.metadata['files'][self.mzml_file]['summary'][fragmentations]['intensity of z=2 phosphoric_acid'] = 'n/a'
                    self.metadata['files'][self.mzml_file]['summary'][fragmentations]['intensity of z=2 water_loss'] = 'n/a'
                    self.metadata['files'][self.mzml_file]['summary'][fragmentations]['z=2 phosphoric_acid to z=2 water_loss intensity ratio'] = 'n/a'
                    self.metadata['files'][self.mzml_file]['summary'][fragmentations]['absolute difference in delta m/z for z=2 phosphoric_acid_loss and z=2 water_loss'] = 'n/a'
            else:
                
                self.metadata['files'][self.mzml_file]['summary'][fragmentations]['intensity of z=2 phosphoric_acid'] = 'n/a'
                self.metadata['files'][self.mzml_file]['summary'][fragmentations]['intensity of z=2 water_loss'] = 'n/a'
                self.metadata['files'][self.mzml_file]['summary'][fragmentations]['z=2 phosphoric_acid to z=2 water_loss intensity ratio'] = 'n/a'
                self.metadata['files'][self.mzml_file]['summary'][fragmentations]['absolute difference in delta m/z for z=2 phosphoric_acid_loss and z=2 water_loss'] = 'n/a'
                self.metadata['files'][self.mzml_file]['summary'][fragmentations]['absolute difference in delta m/z for z=2 phosphoric_acid_loss and z=2 water_loss'] = 'n/a'

            
            #### Detect Phospho or not
            try:
                if self.metadata['files'][self.mzml_file]['summary'][fragmentations]['z=2 phosphoric_acid to z=2 water_loss intensity ratio'] > 1:
                    if self.metadata['files'][self.mzml_file]['neutral_loss_peaks'][loss_type]['phospho_z2']['peak']['mode_bin']['n_spectra'] >= 50:
                        self.metadata['files'][self.mzml_file]['summary'][fragmentations]['has phospho_spectra'] = True
                    elif self.metadata['files'][self.mzml_file]['neutral_loss_peaks'][loss_type]['phosphoric_acid_z2']['peak']['mode_bin']['n_spectra'] >= 50:
                        self.metadata['files'][self.mzml_file]['summary'][fragmentations]['has phospho_spectra'] = True

                    if ('LR' in fragmentations or 'HR_IT_ETD' in fragmentations) and self.metadata['files'][self.mzml_file]['summary'][fragmentations]['absolute difference in delta m/z for z=2 phosphoric_acid_loss and z=2 water_loss'] > 0.3:
                        self.metadata['files'][self.mzml_file]['summary'][fragmentations]['has phospho_spectra'] = False
                    elif 'HR' in fragmentations and self.metadata['files'][self.mzml_file]['summary'][fragmentations]['absolute difference in delta m/z for z=2 phosphoric_acid_loss and z=2 water_loss'] > 0.006:
                        self.metadata['files'][self.mzml_file]['summary'][fragmentations]['has phospho_spectra'] = False
                else:
                    self.metadata['files'][self.mzml_file]['summary'][fragmentations]['has phospho_spectra'] = False
            except:
                self.metadata['files'][self.mzml_file]['summary'][fragmentations]['has phospho_spectra'] = False


            #### Make the call for TMT or iTRAQ
            
            ##Add up all the TMT and iTRAQ values
            for types in results['labeling']['scores']:
                if 'TMT' in types:
                    TMT_all += results['labeling']['scores'][types]
                elif 'iTRAQ' in types:
                    iTRAQ_all += results['labeling']['scores'][types]


            if TMT_all > iTRAQ_all:
                if TMT_all > 2:
                    if results['labeling']['scores']['TMTpro'] > 0.7:
                        results['labeling']['call'] = 'TMTpro'
                    elif results['labeling']['scores']['TMT10'] > 0.7:
                        results['labeling']['call'] = 'TMT10plex'
                    elif results['labeling']['scores']['TMT6'] > 0.7:
                        results['labeling']['call'] = 'TMT6plex'
                    else:
                        results['labeling']['call'] = 'TMT'
                elif results['labeling']['scores']['TMT'] < 1:
                    results['labeling']['call'] = 'none'
                else:
                    results['labeling']['call'] = 'ambiguous'
            else:
                if iTRAQ_all > 0.9:
                    if results['labeling']['scores']['iTRAQ4'] > results['labeling']['scores']['iTRAQ8']:
                        results['labeling']['call'] = 'iTRAQ4plex'
                    else:
                        results['labeling']['call'] = 'iTRAQ8plex'
                elif results['labeling']['scores']['iTRAQ'] < 0.5:
                    results['labeling']['call'] = 'none'
                else:
                    results['labeling']['call'] = 'ambiguous'

            

        self.metadata['files'][self.mzml_file]['summary']['precursor stats'] = self.precursor_stats
        
        full_results = {"fragmentation type": None, "call": None, "fragmentation tolerance": None, "has water_loss":None, "has phospho_spectra": None, "total intensity of z=2 water_loss": 0, "total intensity of z=2 phosphoric_acid": 0, "total z=2 phosphoric_acid to z=2 water_loss intensity ratio":0, "recommended precursor tolerance (ppm)": None}
        self.metadata['files'][self.mzml_file]['summary']['combined summary'] = full_results
        file_summary = self.metadata['files'][self.mzml_file]['summary']
        for keys in file_summary:
            if keys in fragmentation_types:
                try:
                    if full_results["fragmentation type"] is None:
                        full_results["fragmentation type"] = file_summary[keys]["fragmentation_type"]

                    elif full_results["fragmentation type"] != file_summary[keys]["fragmentation_type"]:
                        full_results["fragmentation type"] = "multiple"
                except: 
                    pass
                try:
                    if full_results['call'] is None or full_results['call'] == 'none':
                        full_results["call"] = file_summary[keys]['labeling']["call"]
                    
                    elif file_summary[keys]['labeling']["call"] == 'none':
                        pass

                    elif full_results['call'] != file_summary[keys]['labeling']["call"] and "LR" in keys:
                        pass

                    elif full_results['call'] != file_summary[keys]['labeling']["call"]:
                        full_results['call'] = "multiple"
                except:
                    pass
                try:
                    if full_results['fragmentation tolerance'] is None:
                        full_results['fragmentation tolerance'] = file_summary[keys]['tolerance']
                    elif 'no' in file_summary[keys]['tolerance']:
                        pass
                    elif full_results['fragmentation tolerance'] != file_summary[keys]['tolerance']:
                        full_results['fragmentation tolerance'] = "multiple"

                except:
                    pass

                try:
                    if full_results['has water_loss'] is None:
                        full_results['has water_loss'] = file_summary[keys]['has water_loss']
                    elif file_summary[keys]['has water_loss']:
                        full_results['has water_loss'] = file_summary[keys]['has water_loss']
                except:
                    pass

                try:   
                    if full_results['has phospho_spectra'] is None:
                          full_results['has phospho_spectra'] = file_summary[keys]['has phospho_spectra']
                    elif not full_results['has phospho_spectra'] and file_summary[keys]['has phospho_spectra']:
                        full_results['has phospho_spectra'] = file_summary[keys]['has phospho_spectra']
                    elif full_results['has phospho_spectra']:
                        pass
                    else:
                        full_results['has phospho_spectra'] = False 
                except:
                    pass
                    
                try:
                    
                    full_results['total intensity of z=2 phosphoric_acid'] += file_summary[keys]['intensity of z=2 phosphoric_acid']
                    
                    full_results['total intensity of z=2 water_loss'] += file_summary[keys]['intensity of z=2 water_loss']
                except:
                    pass

                try:
                    full_results['recommended precursor tolerance (ppm)'] = self.precursor_stats["precursor tolerance"]['fit_ppm']['recommended precursor tolerance (ppm)']
                except:
                    pass
        
        full_results['total z=2 phosphoric_acid to z=2 water_loss intensity ratio'] = full_results['total intensity of z=2 phosphoric_acid'] / (full_results['total intensity of z=2 water_loss'] + epsilon)
        

        

    ####################################################################################################
    #### Log an event
    def log_event(self, status, code, message):

        category = 'UNKNOWN'
        if status == 'WARNING':
            category = 'warnings'
        elif status == 'ERROR':
            category = 'errors'
        else:
            eprint(f"FATAL ERROR: Unrecognized event status '{status}'")
            eprint(sys.exc_info())
            sys.exit()

        #### Record the event
        full_message = f"{status}: [{code}]: {message}"
        self.metadata['problems'][category]['count'] += 1
        self.metadata['problems'][category]['list'].append(full_message)
        if code not in self.metadata['problems'][category]['codes']:
            self.metadata['problems'][category]['codes'][code] = 1
        else:
            self.metadata['problems'][category]['codes'][code] += 1

        #### If this is an error, also update the overall state
        if status == 'ERROR':
            self.metadata['state']['status'] = status
            self.metadata['state']['code'] = code
            self.metadata['state']['message'] = message

        #### In verbose mode, print out warnings and errors
        if self.verbose >= 1:
            eprint(full_message)
        else:
            if status == 'ERROR':
                eprint(full_message)

####################################################################################################
#### Gaussian function used for curve fitting
def gaussian_function(x,a,x0,sigma, y_offset):
    return a*numpy.exp(-(x-x0)**2/(2*sigma**2)) + y_offset


####################################################################################################
#### For command-line usage
def main():

    argparser = argparse.ArgumentParser(description='Creates an index for an MSP spectral library file')
    argparser.add_argument('--use_cache', action='store', help='Set to true to use the cached file instead of regenerating')
    argparser.add_argument('--refresh', action='count', default=0, help='If set, existing metadata for a file will be overwritten rather than skipping the file')
    argparser.add_argument('--verbose', action='count', help='If set, print more information about ongoing processing' )
    argparser.add_argument('--version', action='version', version='%(prog)s 0.5')
    argparser.add_argument('files', type=str, nargs='+', help='Filenames of one or more mzML files to read')
    params = argparser.parse_args()

    #### Set verbose
    verbose = params.verbose
    if verbose is None:
        verbose = 1

    #### Loop over all the files to ensure that they are really there before starting work
    for file in params.files:
        if not os.path.isfile(file):
            eprint(f"ERROR: File '{file}' not found or not a file")
            return

    #### Load the current metadata file in order to update it
    study = MetadataHandler(verbose=verbose)
    study.read_or_create()

    #### Loop over all the files processing them
    for file in params.files:

        #### If this file is already in the metadata store, purge it and generate new
        if file in study.metadata['files']:
            if params.refresh is not None and params.refresh > 0:
                study.metadata['files'][file] = None
            else:
                if verbose >= 1:
                    eprint(f"INFO: Already have results for '{file}'. Skipping..")
                continue

        #### Assess the mzML file
        assessor = MzMLAssessor(file, metadata=study.metadata, verbose=verbose)
        assessor.read_header()
        if not params.use_cache:
            assessor.read_spectra()
        assessor.assess_lowend_composite_spectra()
        assessor.assess_neutral_loss_composite_spectra()

        if assessor.metadata['state']['status'] != 'ERROR' or 'multiple fragmentation types' in assessor.metadata['state']['message']:
            assessor.assess_ROIs()
        
    #### Infer parameters based on the latest data
    if assessor.metadata['state']['status'] != 'ERROR' or 'multiple fragmentation types' in assessor.metadata['state']['message']:
        study.infer_search_criteria()

    #### Write out our state of mind
    study.store()


#### For command line usage
if __name__ == "__main__":
    main()

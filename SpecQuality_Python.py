#!/usr/bin/env python3
import sys
import re
import time
import math
import os

#version 0.7
# Amino acid masses
AMINO_ACID_MASSES = {
    'A': 71.03711, 'C': 103.00919, 'D': 115.02694, 'E': 129.04259,
    'F': 147.06841, 'G': 57.02146, 'H': 137.05891, 'I': 113.08406,
    'K': 128.09496, 'L': 113.08406, 'M': 131.04049, 'N': 114.04293,
    'P': 97.05276, 'Q': 128.05858, 'R': 156.10111, 'S': 87.03203,
    'T': 101.04768, 'V': 99.06841, 'W': 186.07931, 'Y': 163.06333,
}

def custom_log10(n):
    """
    Custom log10 to handle non-positive inputs, where subsequent calculations often use fallbacks like 0.01
    """
    if isinstance(n, (int, float)) and n > 0:
        return math.log10(n)
    return -float('inf') # Represents problematic log, handled later

def ppm_to_da(ppm, mass):
    """Converts PPM tolerance to Dalton, ensuring mass is positive."""
    if mass <= 0: # Avoid issues with zero or negative mass in calculation
        return (ppm * 1e-6) # Effectively a very small Da value if mass was meant to be ~1
    return (ppm * mass) / 1_000_000

def calculate_spectral_quality_features(peaks_list, parent_mh_val, precursor_intensity_val, charge_z_val, 
                                        tolerance, tolerance_unit, noise_thresh_percent):
    try:
        original_prec_intensity = float(precursor_intensity_val)
    except ValueError:
        original_prec_intensity = 1.0 # Default

    # precursor intensity is immediately log10 transformed.
    log10_transformed_prec_intensity = custom_log10(original_prec_intensity)

    num_peaks = len(peaks_list)
    if num_peaks == 0: # Should be caught by peak_cutoff, but defensive
        # Return default/neutral values for an empty spectrum
        intensity_feat = log10_transformed_prec_intensity if math.isfinite(log10_transformed_prec_intensity) else 0.01
        return {
            'peak_count': 0.0, 'signal_peaks': 0.0, 'snr': 0.01,
            'intensity': intensity_feat, 'peak_density': 0.01,
            'good_diff_fraction': 0.01, 'complements_fraction': 0.01,
            'isotope_peaks': 0.01, 'neutral_loss_peaks': 0.01,
            'ari': 0.01, 'sqs': None
        }

    # Convert peaks to float [ (mz, intensity), ... ]
    try:
        float_peaks_list = [(float(p[0]), float(p[1])) for p in peaks_list]
    except (ValueError, TypeError) as e:
         # Fallback if peak conversion fails for some reason
        intensity_feat = log10_transformed_prec_intensity if math.isfinite(log10_transformed_prec_intensity) else 0.01
        print(f"Warning: Error converting peak data to float: {e}. Returning default features.")
        return {
            'peak_count': float(num_peaks), 'signal_peaks': 0.0, 'snr': 0.01,
            'intensity': intensity_feat, 'peak_density': 0.01,
            'good_diff_fraction': 0.01, 'complements_fraction': 0.01,
            'isotope_peaks': 0.01, 'neutral_loss_peaks': 0.01,
            'ari': 0.01, 'sqs': None
        }


    max_ms2_intensity = 0.0
    for _, intensity_val in float_peaks_list:
        if intensity_val > max_ms2_intensity:
            max_ms2_intensity = intensity_val

    if max_ms2_intensity == 0: # All peaks have zero intensity
        intensity_feat = log10_transformed_prec_intensity if math.isfinite(log10_transformed_prec_intensity) else 0.01
        peak_density_val = 0.01
        if num_peaks > 0:
            min_mz_spec = min(p[0] for p in float_peaks_list)
            max_mz_spec = max(p[0] for p in float_peaks_list)
            mz_rng = max_mz_spec - min_mz_spec
            if mz_rng > 0 : peak_density_val = (num_peaks * 100.0) / mz_rng
            else: peak_density_val = 0.01 if num_peaks >0 else 0.0
            if peak_density_val <=0 and num_peaks > 0: peak_density_val = 0.01


        return {
            'peak_count': float(num_peaks), 'signal_peaks': 0.0, 'snr': 0.01,
            'intensity': intensity_feat, 'peak_density': peak_density_val,
            'good_diff_fraction': 0.01, 'complements_fraction': 0.01,
            'isotope_peaks': 0.01, 'neutral_loss_peaks': 0.01,
            'ari': 0.01, 'sqs': None
        }

    num_signal_peaks = 0
    sum_signal_intensity = 0.0
    sum_noise_intensity = 1.0  # Init to 1 to prevent div by zero

    min_mz_spec = float('inf')
    max_mz_spec = float('-inf')

    for mz_val, intensity_val in float_peaks_list:
        normalized_ms2_intensity = (intensity_val * 100.0) / max_ms2_intensity
        
        if normalized_ms2_intensity >= noise_thresh_percent:
            num_signal_peaks += 1
            sum_signal_intensity += intensity_val
        else:
            sum_noise_intensity += intensity_val

        if mz_val < min_mz_spec: min_mz_spec = mz_val
        if mz_val > max_mz_spec: max_mz_spec = mz_val
    
    snr_ratio = sum_signal_intensity / sum_noise_intensity if sum_noise_intensity > 0 else sum_signal_intensity / 0.01
    snr_val = custom_log10(snr_ratio)
    if not math.isfinite(snr_val) or snr_val <= 0:
        snr_val = 0.01
        
    mz_range = max_mz_spec - min_mz_spec
    peak_density_val = (num_peaks * 100.0) / mz_range if mz_range > 0 else 0.01
    if peak_density_val <= 0: peak_density_val = 0.01
    
    sum_good_diff_intensity = 0.0
    sum_total_diff_intensity = 0.0
    sum_complement_intensity = 0.0

    for i in range(num_peaks - 1):
        mz_i, intensity_i = float_peaks_list[i]
        for j in range(i + 1, num_peaks):
            mz_j, intensity_j = float_peaks_list[j]
            
            mass_diff_observed = abs(mz_i - mz_j)
            is_good_diff = False

            if 56 <= mass_diff_observed <= 187: # AA mass range
                sum_total_diff_intensity += intensity_i + intensity_j
                
                for aa_mass_theoretical in AMINO_ACID_MASSES.values():
                    current_tol_da = 0.0
                    if tolerance_unit.lower() == 'da':
                        current_tol_da = tolerance
                    elif tolerance_unit.lower() == 'ppm':
                        current_tol_da = ppm_to_da(tolerance, mass_diff_observed)
                    
                    if abs(mass_diff_observed - aa_mass_theoretical) <= current_tol_da:
                        is_good_diff = True
                        break
            
            if is_good_diff:
                sum_good_diff_intensity += intensity_i + intensity_j

            # Complements
            mass_sum_observed = mz_i + mz_j
            # parent_mh_val is $mz*$z from main script
            diff_from_parent_mh = abs(float(parent_mh_val) - mass_sum_observed) 
            
            current_tol_da_comp = 0.0
            if tolerance_unit.lower() == 'da':
                current_tol_da_comp = tolerance
            elif tolerance_unit.lower() == 'ppm':
                # It uses ppm_to_da(tol, diff_from_parent_mh), which means for non-zero
                # diff_from_parent_mh, a match only occurs if tol >= 1,000,000 ppm.
                # Effectively, it requires diff_from_parent_mh to be zero for a PPM match.
                if diff_from_parent_mh == 0: # Perfect match case
                    current_tol_da_comp = 0.0 
                else:
                    current_tol_da_comp = ppm_to_da(tolerance, abs(diff_from_parent_mh))

            if diff_from_parent_mh <= current_tol_da_comp:
                sum_complement_intensity += intensity_i + intensity_j

    good_diff_fraction_val = (100.0 * sum_good_diff_intensity) / sum_total_diff_intensity if sum_total_diff_intensity > 0 else 0.01
    if good_diff_fraction_val <= 0: good_diff_fraction_val = 0.01

    # Complements Fraction calculation
    # $complements_fraction = log10($complements/$intensity); where $intensity is log10_transformed_prec_intensity
    complements_val_for_log = sum_complement_intensity if sum_complement_intensity > 0 else 0.01
    
    complements_fraction_val = 0.01 # Default
    if math.isfinite(log10_transformed_prec_intensity) and log10_transformed_prec_intensity != 0:
        ratio_for_log = complements_val_for_log / log10_transformed_prec_intensity
        complements_fraction_val = custom_log10(ratio_for_log)
    
    if not math.isfinite(complements_fraction_val) or complements_fraction_val == 0 : # Check for 0 as well
         complements_fraction_val = 0.01 # Fallback for -Inf, NaN, or 0

    # Isotope and Neutral Loss Peaks
    # Peaks should be sorted by m/z for this logic to be correct.
    sorted_float_peaks_list = sorted(float_peaks_list, key=lambda p_item: p_item[0])
    num_isotope_peaks = 0
    num_neutral_loss_peaks = 0
    
    for i in range(len(sorted_float_peaks_list) - 1):
        mz1, _ = sorted_float_peaks_list[i]
        mz2, _ = sorted_float_peaks_list[i+1]

        isotope_tol_da_curr = tolerance if tolerance_unit.lower() == 'da' else ppm_to_da(tolerance, mz1)
        if abs((mz2 - mz1) - 1.0) < isotope_tol_da_curr:
            num_isotope_peaks += 1
        
        neutral_loss_masses_options = [17.0, 18.0, 28.0] # NH3, H2O, CO
        for loss_mass in neutral_loss_masses_options:
            nl_tol_da_curr = tolerance if tolerance_unit.lower() == 'da' else ppm_to_da(tolerance, loss_mass)
            if abs((mz2 - mz1) - loss_mass) < nl_tol_da_curr:
                num_neutral_loss_peaks += 1 # A pair can match multiple losses

    isotope_peaks_norm_val = (100.0 * num_isotope_peaks / num_peaks) if num_peaks > 0 else 0.01
    if isotope_peaks_norm_val <= 0: isotope_peaks_norm_val = 0.01

    neutral_loss_peaks_norm_val = (100.0 * num_neutral_loss_peaks / num_peaks) if num_peaks > 0 else 0.01
    if neutral_loss_peaks_norm_val <= 0: neutral_loss_peaks_norm_val = 0.01

    # Average of Relative Peak Intensities (ARI)
    # Perl: $ari = $intensity * $signal_peaks / $peak_count; (where $intensity is log10_transformed_prec_intensity)
    ari_val = 0.01 # Default
    if num_peaks > 0 and math.isfinite(log10_transformed_prec_intensity):
        ari_val = log10_transformed_prec_intensity * num_signal_peaks / num_peaks
    if not math.isfinite(ari_val) or ari_val == 0: ari_val = 0.01

    return {
        'peak_count': float(num_peaks),
        'signal_peaks': float(num_signal_peaks),
        'snr': float(snr_val),
        'intensity': float(log10_transformed_prec_intensity if math.isfinite(log10_transformed_prec_intensity) else 0.01),
        'peak_density': float(peak_density_val),
        'good_diff_fraction': float(good_diff_fraction_val),
        'complements_fraction': float(complements_fraction_val),
        'isotope_peaks': float(isotope_peaks_norm_val),
        'neutral_loss_peaks': float(neutral_loss_peaks_norm_val),
        'ari': float(ari_val),
        'sqs': None
    }

def fix_complements_fraction_values(spectra_data_dict):
    # Modifies 'complements_fraction' in spectra_data_dict in-place.
    if not spectra_data_dict: return

    min_val = float('inf')
    for scan_id in spectra_data_dict:
        features = spectra_data_dict[scan_id]
        val = features.get('complements_fraction')
        if isinstance(val, (int, float)) and math.isfinite(val):
            if val < min_val: min_val = val
    
    if math.isfinite(min_val) and min_val < 0:
        for scan_id in spectra_data_dict:
            features = spectra_data_dict[scan_id]
            current_cf = features.get('complements_fraction')
            if isinstance(current_cf, (int,float)) and math.isfinite(current_cf):
                 features['complements_fraction'] += -min_val + 0.01
            else: # Was NaN, -Inf, or not present; set to a default positive value
                 features['complements_fraction'] = abs(min_val) + 0.01 if math.isfinite(min_val) else 0.01


def calculate_sqs_for_all_spectra(original_spectra_data):
    # Returns a new dictionary with SQS values calculated.
    processed_spectra_data = {}
    features_for_sqs = [
        'peak_count', 'signal_peaks', 'snr', 'intensity', 
        'peak_density', 'good_diff_fraction', 'complements_fraction', 
        'isotope_peaks', 'neutral_loss_peaks', 'ari'
    ]

    for scan_id, features in original_spectra_data.items():
        new_features = features.copy() # Work on a copy
        product_val = 1.0
        num_valid_features = 0
        for feat_name in features_for_sqs:
            val = features.get(feat_name)
            # Ensure value is positive and finite for geometric mean product
            if not (isinstance(val, (int, float)) and math.isfinite(val) and val > 0):
                val = 0.01 # Default small positive if problematic
            product_val *= val
            num_valid_features +=1
        
        geometric_mean_val = product_val ** (1.0 / num_valid_features) if num_valid_features > 0 else 0.0
        
        new_features['sqs'] = float(geometric_mean_val if math.isfinite(geometric_mean_val) else 0.01)
        processed_spectra_data[scan_id] = new_features
        
    return processed_spectra_data

# It also expects a 'minmax_ref' which is not defined/populated in the main script.
# def range_normalize(original_spectra_data, minmax_ref_dict): ...

def main():
    if not (2 <= len(sys.argv) <= 3):
        print(f"Usage: python {os.path.basename(sys.argv[0])} <mgf_file> [msgf_plus_output_file]", file=sys.stderr)
        sys.exit(1)

    mgf_filepath = sys.argv[1]
    msgf_filepath = sys.argv[2] if len(sys.argv) > 2 else None

    # Parameters
    param_tolerance = 0.5
    param_tolerance_unit = 'Da'
    param_noise_threshold = 5.0  # % normalized intensity
    param_peak_cutoff = 7

    script_start_time = time.time()

    # Generate output filenames
    base_mgf_name = mgf_filepath
    if mgf_filepath.lower().endswith('.mgf'):
        base_mgf_name = mgf_filepath[:-4]
    
    sqs_output_filepath = f"{base_mgf_name}_SQS.tsv"
    psm_sqs_output_filepath = f"{base_mgf_name}_PSM_SQS.tsv"

    all_spectra_features = {} 
    
    # MGF parsing state variables
    current_spectrum_peaks = []
    current_title, current_rt_seconds, current_pepmass_mz_str, current_charge_str, current_pepmass_intensity_str = \
        None, None, None, None, None
    default_precursor_intensity = 1.0
    mgf_spectra_count = 0

    print(f"Processing MGF file: {mgf_filepath}")
    try:
        with open(mgf_filepath, 'r') as mgf_file:
            for line_num, line_content in enumerate(mgf_file, 1):
                line = line_content.strip()
                if not line: continue # Skip empty lines

                if line.startswith('BEGIN IONS'):
                    current_spectrum_peaks = []
                    current_title, current_rt_seconds, current_pepmass_mz_str, \
                    current_charge_str, current_pepmass_intensity_str = '', '', '', '', ''
                elif line.startswith('TITLE='):
                    current_title = line[len('TITLE='):]
                elif line.startswith('RTINSECONDS='):
                    current_rt_seconds = line[len('RTINSECONDS='):]
                elif line.startswith('PEPMASS='):
                    pepmass_data = line[len('PEPMASS='):]
                    parts = pepmass_data.split()
                    current_pepmass_mz_str = parts[0]
                    if len(parts) > 1:
                        current_pepmass_intensity_str = parts[1]
                        try: # Update default intensity if provided and valid
                            default_precursor_intensity = float(current_pepmass_intensity_str)
                        except ValueError: # If intensity is not a valid float, keep old default
                            pass 
                    else:
                        current_pepmass_intensity_str = str(default_precursor_intensity)
                elif line.startswith('CHARGE='):
                    current_charge_str = line[len('CHARGE='):].replace('+', '')
                elif re.match(r"^\d+\.?\d*\s+\d+\.?\d*", line): # Peak data
                    try:
                        mass_str, intensity_str = line.split()
                        current_spectrum_peaks.append((mass_str, intensity_str)) # Keep as str for now
                    except ValueError:
                        print(f"Warning line {line_num}: Could not parse peak data '{line}' for '{current_title}'. Skipping peak.", file=sys.stderr)
                elif line.startswith('END IONS'):
                    mgf_spectra_count += 1
                    if not all([current_title, current_pepmass_mz_str, current_charge_str]):
                        print(f"Warning line {line_num}: Missing essential info (TITLE, PEPMASS, or CHARGE) for spectrum '{current_title}'. Skipping.", file=sys.stderr)
                        continue
                    
                    if len(current_spectrum_peaks) < param_peak_cutoff:
                        continue

                    try:
                        # Convert necessary precursor info to float for calculations
                        pepmass_mz_float = float(current_pepmass_mz_str)
                        charge_float = float(current_charge_str) # Can be non-integer
                        pepmass_intensity_float = float(current_pepmass_intensity_str)
                        rt_float = float(current_rt_seconds) if current_rt_seconds else 0.0
                    except ValueError as e:
                        print(f"Warning line {line_num}: Cannot convert precursor info to numeric for '{current_title}': {e}. Skipping.", file=sys.stderr)
                        continue
                    
                    # $precursorMH = ($mz*$z); This is not [M+H]+, but used as 'parent_MH' in calculations
                    calculated_parent_mh = pepmass_mz_float * charge_float

                    features_dict = calculate_spectral_quality_features(
                        current_spectrum_peaks, calculated_parent_mh, pepmass_intensity_float,
                        charge_float, param_tolerance, param_tolerance_unit, param_noise_threshold
                    )
                    
                    features_dict['COL'] = [current_pepmass_mz_str, current_charge_str, current_rt_seconds or '0']
                    all_spectra_features[current_title] = features_dict
    except FileNotFoundError:
        print(f"Error: MGF file not found: {mgf_filepath}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred during MGF processing: {e}", file=sys.stderr)
        sys.exit(1)

    mgf_processing_time = time.time()
    print(f"MGF file processed in: {mgf_processing_time - script_start_time:.2f} secs ({(mgf_processing_time - script_start_time)/60:.2f} mins)")

    fix_complements_fraction_values(all_spectra_features)
    
    print("Calculating SQS...")
    all_spectra_features_with_sqs = calculate_sqs_for_all_spectra(all_spectra_features)
    print(f"SQS calculated for {len(all_spectra_features_with_sqs)} spectra")

    sqs_output_header = "scan\tmz\tz\tRT\tPeakCount\tSignalPeaksCount\tSNR\tIntensity\tPeakDensity\tGoodDiffFraction\tComplementsFraction\tIsotopePeaks\tNeutralLossPeaks\tAverageRelativeIntensity\tSQS_GM\n"
    try:
        with open(sqs_output_filepath, 'w') as sqs_out_file:
            sqs_out_file.write(sqs_output_header)
            for scan_id, features in all_spectra_features_with_sqs.items():
                col_data = features.get('COL', ['NA', 'NA', 'NA'])
                row = [
                    scan_id, col_data[0], col_data[1], col_data[2],
                    f"{features.get('peak_count', 0.0):.4f}", f"{features.get('signal_peaks', 0.0):.4f}",
                    f"{features.get('snr', 0.0):.4f}", f"{features.get('intensity', 0.0):.4f}",
                    f"{features.get('peak_density', 0.0):.4f}", f"{features.get('good_diff_fraction', 0.0):.4f}",
                    f"{features.get('complements_fraction', 0.0):.4f}", f"{features.get('isotope_peaks', 0.0):.4f}",
                    f"{features.get('neutral_loss_peaks', 0.0):.4f}", f"{features.get('ari', 0.0):.4f}",
                    f"{features.get('sqs', 0.0):.4f}"
                ]
                sqs_out_file.write("\t".join(row) + "\n")
    except IOError as e:
        print(f"Error writing SQS output file {sqs_output_filepath}: {e}", file=sys.stderr)
        sys.exit(1)
    
    sqs_writing_time = time.time()
    print(f"SQS calculated and Features printed in: {sqs_writing_time - mgf_processing_time:.2f} secs ({(sqs_writing_time - mgf_processing_time)/60:.2f} mins)")

    if msgf_filepath:
        print(f"Processing MSGF+ file: {msgf_filepath}")
        if not os.path.exists(msgf_filepath):
            print(f"Warning: MSGF+ file '{msgf_filepath}' not found. Skipping merge step.", file=sys.stderr)
        else:
            msgf_header_written = False
            msgf_matched_count = 0
            msgf_unmatched_count = 0
            msgf_total_scans = 0
            try:
                with open(msgf_filepath, 'r') as msgf_in_file, open(psm_sqs_output_filepath, 'w') as psm_out_file:
                    for line in msgf_in_file:
                        original_line = line.rstrip('\n')
                        if not msgf_header_written:
                            psm_out_file.write(f"{original_line}\tPeakCount\tSignalPeaksCount\tSNR\tIntensity\tPeakDensity\tGoodDiffFraction\tComplementsFraction\tIsotopePeaks\tNeutralLossPeaks\tAverageRelativeIntensity\tSQS\tLabel\n")
                            msgf_header_written = True
                            continue
                        
                        msgf_total_scans += 1
                        parts = original_line.split('\t')
                        # Perl: $scan=$arr[3]; $qval=$arr[16]; (0-indexed from split)
                        # Adjust indices if your MSGF+ TSV format differs.
                        if len(parts) < 17: 
                            print(f"Warning: MSGF+ line too short: '{original_line[:50]}...'. Skipping.", file=sys.stderr)
                            msgf_unmatched_count +=1
                            continue
                        
                        scan_title_from_msgf = parts[3]
                        q_value_str = parts[16]
                        
                        try:
                            q_value_float = float(q_value_str)
                        except ValueError:
                            q_value_float = float('inf') # Treat unparseable Q-value as poor

                        if scan_title_from_msgf in all_spectra_features_with_sqs:
                            msgf_matched_count += 1
                            sqs_features = all_spectra_features_with_sqs[scan_title_from_msgf]
                            
                            label = 'NA'
                            if q_value_float <= 0.01: label = 'Excellent'
                            elif q_value_float <= 0.05: label = 'Good'
                            elif q_value_float <= 0.1: label = 'Average'
                            else: label = 'Poor'

                            features_to_append = [
                                f"{sqs_features.get('peak_count', 0.0):.4f}", f"{sqs_features.get('signal_peaks', 0.0):.4f}",
                                f"{sqs_features.get('snr', 0.0):.4f}", f"{sqs_features.get('intensity', 0.0):.4f}",
                                f"{sqs_features.get('peak_density', 0.0):.4f}", f"{sqs_features.get('good_diff_fraction', 0.0):.4f}",
                                f"{sqs_features.get('complements_fraction', 0.0):.4f}", f"{sqs_features.get('isotope_peaks', 0.0):.4f}",
                                f"{sqs_features.get('neutral_loss_peaks', 0.0):.4f}", f"{sqs_features.get('ari', 0.0):.4f}",
                                f"{sqs_features.get('sqs', 0.0):.4f}", label
                            ]
                            psm_out_file.write(f"{original_line}\t" + "\t".join(features_to_append) + "\n")
                            # Perl deletes scan from %sp; not strictly needed if MSGF+ scan titles are unique
                        else:
                            msgf_unmatched_count += 1
                            # Optionally write line with NA for SQS features if all MSGF+ lines must be in output
                            # psm_out_file.write(f"{original_line}\t" + "\t".join(['NA'] * 12) + "\n")


                print(f"MSGF+ scans (Total): {msgf_total_scans}, Matched: {msgf_matched_count}, Unmatched: {msgf_unmatched_count}")
            except IOError as e:
                print(f"Error processing MSGF+ file or writing PSM SQS output: {e}", file=sys.stderr)
            except Exception as e:
                print(f"An unexpected error occurred during MSGF+ processing: {e}", file=sys.stderr)

    script_end_time = time.time()
    total_duration = script_end_time - script_start_time
    print(f"Total time taken: {total_duration:.2f} secs ({total_duration/60:.2f} mins)")

if __name__ == '__main__':
    main()

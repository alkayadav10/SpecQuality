# SpecQuality
SpecQuality is a tool for quality assessment of MS/MS spectra
## Spectral Quality Score (SQS) Calculator

**Version:** 1.0

**Language:** Python / Perl

### 1. Documentation

#### 1.1. Purpose

This tool processes tandem mass spectrometry (MS/MS) data from MGF files to assess the quality of each spectrum. It calculates a set of ten distinct spectral features and then combines them into a **Spectral Quality Score (SQS)**, which are defined based on how they are calculated:-
  * the geometric mean of 5 features (**SQS5_gm**)
  * the geometric mean of 10 features (**SQS5_XGB_clf**)
  * 5 feature based XGBoost classifier trained SQS (**SQS5_XGB_clf**)
  * 10 feature based XGBoost classifier trained SQS (**SQS10_XGB_clf**)
  * 5 feature based XGBoost regression trained SQS (**SQS5_XGB_reg**)
  * 5 feature based XGBoost regression trained SQS (**SQS5_XGB_reg**)

**Recommendation**:
  * **SQS10_gm**
  * **SQS5_gm** are recommended, followed by
  * **SQS10_XGB_clf**
  * **SQS5_XGB_clf** 


Optionally, if an MSGF+ (a common peptide identification search engine) output file (typically in TSV format) is provided, the script will merge the calculated SQS and its component features with the Peptide-Spectrum Matches (PSMs) from the MSGF+ results. This allows for correlating spectral quality with identification confidence (e.g., Q-value).

Also check project Wiki [https://github.com/alkayadav10/SpecQuality/wiki/Documentation] for details

# SpecQuality: Spectral Quality Prediction Tool

A comprehensive tool for predicting spectral quality in mass spectrometry data using both geometric mean methods and machine learning approaches.

## Overview

SpecQuality calculates spectral quality scores (SQS) for mass spectrometry data through two complementary approaches:

1. **Traditional Method**: Geometric mean of 10 spectral features
2. **Machine Learning Method**: XGBoost-based prediction using pre-trained models

The tool supports both **MGF** and **mzML** file formats and provides comprehensive feature extraction and quality assessment.

## Quick Start

### Basic Feature Extraction (Geometric Mean Only)

```bash
# For MGF files:
python SpecQuality_Calculate_Features.py input_file.mgf

# For mzML files:
python SpecQuality_Calculate_Features.py input_file.mzml

#### 1.2. Key Features Calculated


For each spectrum, the following features are calculated:

 

1.  **Peak Count (`peak_count`):** Total number of peaks in the MS/MS spectrum.

2.  **Signal Peaks Count (`signal_peaks`):** Number of peaks whose intensity is above a defined noise threshold (relative to the most intense peak in the spectrum).

3.  **Signal-to-Noise Ratio (`snr`):** Log10 of the ratio of total intensity of signal peaks to total intensity of noise peaks.

4.  **Intensity (`intensity`):** Log10 of the precursor ion intensity (if available in MGF, otherwise defaults to TIC based Intensity summed up from MS/MS intensities). 

5.  **Peak Density (`peak_density`):** Number of peaks per 100 m/z range.

6.  **Good Difference Fraction (`good_diff_fraction`):** Percentage of intensity from peak pairs whose m/z difference corresponds to a known amino acid mass within a given tolerance.

7.  **Complements Fraction (`complements_fraction`):** A measure related to the intensity of peak pairs whose m/z sum is close to the precursor m/z multiplied by its charge (a proxy for parent MH). This value is log-transformed and then normalized across all spectra to ensure positivity.

8.  **Isotope Peaks (`isotope_peaks`):** Normalized count of adjacent peaks separated by ~1 Da (indicative of isotopic distributions).

9.  **Neutral Loss Peaks (`neutral_loss_peaks`):** Normalized count of adjacent peaks whose m/z difference corresponds to common neutral losses (e.g., H2O, NH3, CO).

10. **Average Relative Intensity (`ari`):** Product of the log10 precursor intensity and the ratio of signal peaks to total peaks.

**Spectral Quality Score (`sqs`):** The geometric mean of the ten features listed above.



 

#### 1.3. Input Files

 

1.  **Spectra File (Required): (mzML or MGF format)**

    *   A standard HUPO-PSI mzML format, or
    *   Mascot Generic Format file containing MS/MS spectra.

    *   Expected fields per spectrum for MGF:

        *   `BEGIN IONS` / `END IONS`

        *   `TITLE=` (Used as scan identifier)

        *   `PEPMASS=` (Precursor m/z, optionally precursor intensity)

        *   `CHARGE=` (Precursor charge state)

        *   `RTINSECONDS=` (Retention time, optional)

        *   Peak list (m/z and intensity pairs)

2.  **MSGF+ Output File (Optional):**

    *   A tab-separated values (TSV) file, typically an output from the MSGF+ search engine.

    *   The script expects to find the scan title/identifier (MGF compatible) in the 4th column (index 3) and the Q-value in the 17th column (index 16) of this file. These indices might need adjustment if your MSGF+ output format differs.

 

#### 1.4. Output Files

 

The script generates one or two output files in the same directory as the input MGF file:

 

1.  **`<msmsinput_file_base>_SQS.tsv`:**

    *   A tab-separated file containing the calculated spectral features and SQS for every spectrum in the input MGF file that passes the minimum peak cutoff.

    *   Columns: `scan`, `mz`, `z`, `RT`, `PeakCount`, `SignalPeaksCount`, `SNR`, `Intensity`, `PeakDensity`, `GoodDiffFraction`, `ComplementsFraction`, `IsotopePeaks`, `NeutralLossPeaks`, `AverageRelativeIntensity`, `SQS`.

 

2.  **`<msms_file_base>_PSM_SQS.tsv` (Generated if MSGF+ file is provided):**

    *   A tab-separated file that merges the original lines from the MSGF+ output with the corresponding spectral quality features and spectral quality score (SQS).

    *   It appends the quality features/SQS and a `Label` column (categorizing PSM quality based on Q-value: "Excellent", "Good", "Average", "Poor") to each matched MSGF+ entry.

 

#### 1.5. Internal Parameters & Behavior

 

*   **Tolerance for mass matching:** `0.5 Da` (hardcoded, `param_tolerance`, `param_tolerance_unit`). This is used for amino acid mass differences, complement matching, isotope, and neutral loss detection. This should be edited to match the experimental/instrumental tolerance.

*   **Noise Threshold:** `5.0%` (hardcoded, `param_noise_threshold`). Peaks with normalized intensity below this threshold are considered noise.

*   **Minimum Peak Cutoff:** `7` (hardcoded, `param_peak_cutoff`). Spectra with fewer than this many peaks are skipped.

*   **Logarithm Handling:** A custom `custom_log10` function is used, which returns `-inf` for non-positive inputs. Subsequent calculations often use a fallback value (e.g., `0.01`) if a log result is not finite or positive.

*   **Complements Fraction Normalization:** The `complements_fraction` feature undergoes a post-processing step (`fix_complements_fraction_values`) where the minimum negative value across all spectra is identified, and this (plus a small epsilon of `0.01`) is added to all `complements_fraction` values to ensure they are positive before SQS calculation. This is crucial for the geometric mean.

*   **PPM to Dalton Conversion:** `ppm_to_da` function converts a PPM tolerance to Daltons based on a given mass. Note specific handling for complement fraction PPM calculations.

 

### 2. Installation Guide

 

#### 2.1. Requirements

 

*   **Python 3:** This script is written for Python 3 (e.g., Python 3.6 or newer recommended).

*   **Standard Python Libraries:** The script uses only standard Python libraries (`sys`, `re`, `time`, `math`, `os`), so no additional package installations are typically required if you have a standard Python 3 environment.

 

#### 2.2. Setup

 

1.  **Ensure Python 3 is installed:**

    Open a terminal or command prompt and type:

    ```bash

    python3 --version

    ```

    If Python 3 is not installed, download it from [python.org](https://www.python.org/downloads/).

 

2.  **Save the Script:**

    Save the provided Python code to a file. For example, name it `sqs_calculator.py`.

 

3.  **Make the Script Executable (Optional, for Linux/macOS):**

    In your terminal, navigate to the directory where you saved the file and run:

    ```bash

    chmod +x sqs_calculator.py

    ```

    This allows you to run it directly as `./sqs_calculator.py` instead of `python3 sqs_calculator.py`.

 

### 3. Quick Start Guide

 

#### 3.1. Prepare Your Input Files

 

1.  **MGF File:** Have your MGF file ready (e.g., `my_experiment.mgf`).

2.  **MSGF+ File (Optional):** If you want to merge with PSM data, have your MSGF+ output TSV file ready (e.g., `my_experiment_msgfplus.tsv`). Ensure the scan identifier is in the 4th column and Q-value in the 17th.

 

#### 3.2. Running the Script

 

Open your terminal or command prompt and navigate to the directory containing the script and your input files.

 

**Scenario 1: Calculate SQS from MGF only**

 

```bash

python3 sqs_calculator.py my_experiment.mgf

```

Or, if you made it executable (Linux/macOS):

```bash

./sqs_calculator.py my_experiment.mgf

```

 

**Scenario 2: Calculate SQS and merge with MSGF+ results**

 

```bash

python3 sqs_calculator.py my_experiment.mgf my_experiment_msgfplus.tsv

```

Or, if executable:

```bash

./sqs_calculator.py my_experiment.mgf my_experiment_msgfplus.tsv

```

 

#### 3.3. Check the Output

 

After the script finishes processing, you will find new `.tsv` files in the same directory:

 

*   `my_experiment_SQS_V7.tsv`: Contains the SQS and detailed features for each spectrum.

*   `my_experiment_PSM_SQS_V7.tsv` (if MSGF+ file was provided): Contains your MSGF+ data augmented with the SQS features and a quality label.

 

You can open these TSV files with spreadsheet software (like Excel, LibreOffice Calc) or text editors that handle tab-separated data.

 

**Example Output Snippet (`_SQS_V7.tsv`):**

 

```

scan    mz      z       RT      PeakCount       SignalPeaksCount SNR     Intensity       PeakDensity     GoodDiffFraction        ComplementsFraction     IsotopePeaks    NeutralLossPeaks        AverageRelativeIntensity        SQS_GM

ControllerX.scanY.scanY.2       600.1234        2       1800.5  50.0000 25.0000 1.5000  6.0000  5.0000  15.0000 0.5000  10.0000 5.0000  3.0000  3.4567

... (more spectra) ...

```

 

This should give you a good starting point for using the script and understanding its outputs.

---

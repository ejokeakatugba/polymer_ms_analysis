# Polymer MS Analysis

**Python-Based Tools for Polymer and Copolymer Mass Spectrometry Analysis**

## Overview

Polymer mass spectra often contain complex distributions of ions arising from repeating units, oligomeric species, and different adduct formations. Identifying these patterns and extracting meaningful chemical information can be challenging, particularly for copolymers containing more than one repeating unit.

This repository contains Python-based computational tools developed to simplify polymer mass spectrometry data analysis through automated repeating-unit identification, spectral signal analysis, and molecular formula calculations.

The tools were developed as part of research conducted in the Gamez Research Group at Texas Tech University.

## Research Software

### 1. Polymer Repeating-Unit Analysis

**File:** `Polymer_analysis.py`

**Developers:** Dr. Dong Zhang and Ejoke Akatugba

Originally developed in connection with our research on single-particle microplastics using laser-assisted micro-pyrolysis flowing atmospheric-pressure afterglow mass spectrometry imaging (LAMP-FAPA MSI).

**Key Features:**
- Searches mass spectra for characteristic polymer repeating-unit patterns.
- Identifies m/z signals separated by specified repeating-unit masses.
- Counts detected signals within identified polymer series.
- Calculates summed signal intensities associated with the detected series.
- Supports polymer identification and intensity-based comparisons.

### 2. Copolymer Mass Spectrometry Analysis

**File:** `copolymer4v3.py`

**Developer:** Ejoke Akatugba

Independently developed to extend repeating-unit analysis to copolymers containing two chemically distinct repeating units.

The current implementation focuses on polyethylene glycol–polyethylene (PEG–PE) copolymers.

**Key Features:**
- Searches for two different polymer repeating-unit series within the same mass spectrum.
- Identifies PEG (44.026215 Da) and PE (28.031300 Da) repeating-unit patterns.
- Applies a dynamic mass tolerance of 5 ppm.
- Calculates summed intensities of identified polymer series.
- Evaluates relationships between series to identify candidate copolymer patterns.
- Exports organized analysis results to Excel.

### 3. Copolymer Exact Molecular Formula Calculator

**File:** `copolymer_calculator.py`

**Developer:** Ejoke Akatugba

Independently developed to support molecular formula interpretation of copolymer ions using combinations of two repeating units and possible adduct formations.

**Key Features:**
- Calculates candidate molecular formulas containing two distinct polymer repeating units.
- Accounts for combinations of repeating units within copolymer compositions.
- Considers possible adduct formations.
- Calculates theoretical exact masses to support candidate ion assignments.

The calculator is intended to assist molecular formula interpretation; candidate assignments require confirmation using appropriate experimental evidence.

## Scientific Significance

These computational tools address challenges associated with interpreting complex polymer mass spectra, particularly the identification of repeating-unit patterns and relationships between oligomer series.

By reducing repetitive manual searches and organizing related spectral signals, the programs support more systematic polymer characterization.

The progression from single-polymer analysis to dual repeating-unit searches and molecular formula calculations extends these capabilities toward the interpretation of more complex copolymer systems.

## Requirements

- Python 3.11 or later
- NumPy
- SciPy
- Pandas
- OpenPyXL
- XlsxWriter

Install the existing dependency list using `pip install -r requirement.txt`.

For the copolymer program's Excel output, also install `xlsxwriter` using `pip install xlsxwriter`.

Additional requirements for the formula calculator should be checked against its imports.

## Usage

**Polymer repeating-unit analysis**

Run `python Polymer_analysis.py`.

Update the input Excel file path and analysis parameters within the script before execution.

**Copolymer analysis**

Run `python copolymer4v3.py`.

Select an Excel file containing mass-to-charge ratios and signal intensities when prompted.

**Exact molecular formula calculator**

Run `python copolymer_calculator.py`, provided the uploaded script uses this filename and supports direct execution.

## Associated Publication

**Single-Particle Microplastics Analysis via Laser-Assisted MicroPyrolysis Flowing Atmospheric-Pressure Afterglow Mass Spectrometry Imaging**

*Analytical Chemistry* (2026)

https://doi.org/10.1021/acs.analchem.6c02632

## Authors and Contributions

**Dr. Dong Zhang**
- Co-developer of the original polymer repeating-unit analysis program.

**Ejoke Akatugba**
- Co-developer of the original polymer repeating-unit analysis program.
- Independent developer of the copolymer mass spectrometry analysis program.
- Independent developer of the copolymer exact molecular formula calculator.

## Affiliation

Gamez Research Group  
Department of Chemistry and Biochemistry  
Texas Tech University

# ESRA

## Overview

ESRA is an interactive Streamlit application developed to identify peptide subsequences that satisfy a consensus score derived from computational peptide analysis.

Starting from a protein (or peptide sequence), the software:

* Generates all unique subsequences of length 7–8 amino acids.
* Computes physicochemical descriptors using the Python package `peptides`.
* Applies classification rules.
* Scores each peptide according to the number of satisfied rules.
* Produces a ranked report highlighting the most promising candidates.
* Exports results as an Excel file.

---

## Associated Publication

This software accompanies the following publication:

**Davide Sestaioni, Simone Ventisette, Giulia Ciacci, Pasquale Palladino, Andrea Barucci, Maria Minunni, Simona Scarano**

*Antibody-Free SPR Detection of Human Myoglobin in Serum by a Sandwich Configuration of Epitope-Imprinted Nanofilms and Nanoparticles*

ACS Sensors, 2026

DOI: DOI_TO_BE_INSERTED

---

## Installation

### 1. Install Python

Python 3.10 or newer is recommended.

Download Python from:

https://www.python.org/downloads/

---

### 2. Download the repository

Clone the repository:

```bash
git clone https://github.com/GAIA-IFAC-CNR/ESRA.git
cd ESRA
```

Alternatively, download the repository as a ZIP file from GitHub.

---

### 3. Install dependencies

Open a terminal inside the project folder and run:

```bash
pip install -r requirements.txt
```

---

## Launching the Application

From the project directory execute:

```bash
streamlit run ESRA.py
```

A browser window will open automatically.

If it does not open, copy the local URL shown in the terminal and paste it into your browser.

---

## Using the Application

### Step 1

Insert:

* Sequence name
* Amino acid sequence


### Step 2

Click: **Generate Report**

### Step 3

Inspect the generated report.

### Step 4

Download the Excel file using the provided download button.

### Step 5

To stop the execution press ctrl+C inside the terminal.

---

## Output

The generated report contains:

* Peptide identifier
* Peptide sequence
* Consensus score
* Descriptor values
* Failed descriptors
* Visual highlighting of rule violations

Only peptides satisfying at least 5 consensus rules are retained.

---

## Example

An example sequence is available in:

```text
examples/example_sequence.txt
```

An example report is available in:

```text
examples/example_output.xlsx
```

---

## Citation

If you use this software in your work, please cite:

**Davide Sestaioni, Simone Ventisette, Giulia Ciacci, Pasquale Palladino, Andrea Barucci, Maria Minunni, Simona Scarano**

*Antibody-Free SPR Detection of Human Myoglobin in Serum by a Sandwich Configuration of Epitope-Imprinted Nanofilms and Nanoparticles*

ACS Sensors, 2026

DOI: DOI_TO_BE_INSERTED

---

## Contact

Principal investigator: 

**Simona Scarano**

Department of Chemistry “Ugo Schiff’, University of Florence, Via della Lastruccia, 3-13, 50019, Sesto Fiorentino, Italy

Email: simona.scarano@unifi.it


---

## Affiliations

Department of Chemistry “Ugo Schiff’, University of Florence, Via della Lastruccia, 3-13, 50019, Sesto Fiorentino, Italy

 Institute of Applied Physics “Nello Carrara”, National Research Council, Via Madonna del Piano 10, Sesto Fiorentino, 50019, Florence, Italy

Department of Pharmacy, University of Pisa, Via Bonanno, 6, 56126, Pisa, Italy


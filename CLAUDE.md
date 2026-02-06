# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository contains tools for parsing and analyzing PacBio Sequel IIe sequencing run metrics. It extracts run conditions and statistics from XML files and performs statistical analysis of loading metrics (P0, P1, P2) vs yield at the SMRT Cell level.

## Commands

### Parse Sequel Run Data (Python)
```bash
python parse_sequel_runs.py /path/to/sequel/data -o output.csv --excel
```
- Requires pandas and openpyxl (for Excel output)
- Scans directory recursively for `.consensusreadset.xml`, `.subreadset.xml`, and `.sts.xml` files
- Filters out barcode-specific files to provide run-wide metrics only

### Run Statistical Analysis (R)
```bash
Rscript analyze_pacbio_runs.R
```
- Expects `ALL_260206.csv` (or similar) in the working directory
- Outputs plots to `figures/` and CSVs to `results/`
- Required R packages: tidyverse, ggplot2, ggpubr, scales, broom

## Architecture

### Data Pipeline
1. **parse_sequel_runs.py** - Main parser that extracts metrics from PacBio XML files
   - `SequelRunParser` class handles all XML parsing with PacBio-specific namespaces
   - Handles both CCS/HiFi (ConsensusReadSet) and CLR (SubreadSet) run types
   - Extracts: run conditions (instrument, kits, movie length, insert size) and statistics (yield, P0/P1/P2 loading, read lengths, SNR)

2. **analyze_pacbio_runs.R** - Statistical analysis and visualization
   - Correlation analysis between loading metrics and yield
   - Linear regression models comparing P1-only vs multi-predictor models
   - Generates distribution histograms, scatter plots, and diagnostic plots

### Supporting Debug Scripts
- `inspect_xml.py` - Dumps XML structure of `.sts.xml` files
- `debug_sts.py` - Tests STS file parsing
- `check_xml_versions.py` - Identifies different XML schema versions across files

## Key Domain Concepts

- **P0/P1/P2**: ZMW loading metrics (Empty/Single/Multi polymerase loading)
- **SMRT Cell**: Individual sequencing chip; each row in output = one SMRT Cell
- **CCS/HiFi vs CLR**: Two sequencing modes with different file formats
- **Yield (Gb)**: Total sequencing output in gigabases
- **Insert Size**: Target library fragment size
- **total_length**: Sum of HiFi consensus read lengths (from consensusreadset.xml)
- **total_bases**: Raw polymerase sequencing output (from sts.xml SequencingUmy)

## Current Issue: XML TotalLength Parsing

Some runs show abnormally high `total_length` values (100-137 Gb vs typical 20-80 Gb). Symptoms:
- `total_length / num_records` gives 19-24 kb avg read length (should be ~insert_size, 8-11 kb)
- `total_length / total_bases` ratio is ~0.99 (should be ~0.3-0.5 for HiFi)

**Suspected cause**: Parser may be grabbing wrong `<TotalLength>` element. PacBio XML has multiple:
1. `<DataSetMetadata><TotalLength>` - dataset-level total (CORRECT)
2. `<ExternalResources><ExternalResource>...<TotalLength>` - per-BAM file totals

**Fix attempted** (commit 6c17904): Changed XPath from `.//pbds:TotalLength` (anywhere) to `root.find('pbds:DataSetMetadata').find('pbds:TotalLength')` (specific path).

**Next steps** - inspect XML structure of a high-yield run:
```bash
# Example high-yield file to inspect:
XML_FILE="/projects/codon_0000/data/RapunzlSequel2e/r64241e_20240622_072411/4_D01/m64241e_240625_182743.consensusreadset.xml"

# Check all TotalLength elements and their paths:
grep -n "TotalLength" "$XML_FILE"

# Or use Python to inspect structure:
python3 -c "
import xml.etree.ElementTree as ET
tree = ET.parse('$XML_FILE')
root = tree.getroot()
print('Root tag:', root.tag)
for elem in root.iter():
    if 'TotalLength' in elem.tag:
        print(f'  {elem.tag}: {elem.text}')
"
```

**Root cause found**: The XML files reference different BAM types:
- `.hifi_reads.bam` = Q20+ filtered HiFi reads (correct yield, ~30% of raw)
- `.reads.bam` = ALL CCS reads including low-quality (inflated yield, ~99% of raw)

The parser now extracts `bam_file` and `is_hifi_bam` fields. The R script filters:
```r
filter(is_hifi_bam == TRUE)  # Only include runs with .hifi_reads.bam
```

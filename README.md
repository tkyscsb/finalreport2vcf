# FinalReport to VCF Converter

A fast, multithreaded Python tool to convert Illumina GenomeStudio FinalReport (Matrix format) to chromosome-separated VCF.gz files,
 specifically designed for the **Infinium Asian Screening Array (ASA-24v1-0)**.

## Features

- **Strand-aware**: Correctly handles Illumina TOP/BOT strands against reference FASTA.
- **Locus Arbitration**: Safely resolves duplicate probes / multi-allelic calls for the same locus based on empirically verified log
ic.
- **Multithreaded**: Processes autosomes (Chr 1–22) in parallel with low memory footprint (~3 GB for 1,700 samples).
- **Quality Preservation**: Retains GenCall Score (`GCSCORE`) in VCF `FORMAT` fields.

## GenomeStudio FinalReport Export Settings

This script expects a Matrix-format FinalReport generated with the following GenomeStudio settings:

- **Samples**: All Samples
- **Options**:
  - `use the redo with the best: 50th Percentile GC Score`
  - `Remove zeroed SNPs from Report`: Checked
  - `Exclude intensity only SNPs`: Checked
- **Output Format**: Matrix
- **Strand**: Use Top Strand
- **Scores**: Include GenCall score
- **Delimiter**: Comma (CSV)
- **Map Files**: Do not create map files

## Installation

```bash
git clone [https://github.com/](https://github.com/)tkyscs/finalreport2vcf.git
cd finalreport2vcf
pip install -r requirements.txt
```

## Usage
```bash
python finalreport_to_vcf.py \
  --manifest ASA-24v1-0_A1.csv \
  --refseq /path/to/Homo_sapiens_assembly19.fasta \
  --final-report /path/to/FinalReport.csv \
  --out-prefix my_cohort \
  --threads 8
```

<p align="center">
  <a href="https://github.com/PathoGenOmics-Lab/fstigo">
    <img src="https://github.com/PathoGenOmics-Lab/fstigo/blob/main/.github/logs/fstigo.png" height="250" alt="fstigo">
  </a>
</p>

# FSTiGO
FSTiGO is a command-line tool that calculates pairwise Fst values between samples using allele frequencies from VCF files. It uses a multi-threaded approach to speed up the computation of Fst values, making it well-suited for analyzing large genomic datasets.

## Features

Efficient calculation of pairwise Fst values between multiple samples.

Supports input in the form of VCF files or a list of VCF files.

Uses multi-threading with Rayon for parallel computation.

Displays a progress bar for tracking the computation status.


## Usage

To run FSTiGO, use the following command:
```
fstigo -v <VCF_FILES> -r <REFERENCE_FASTA> -o <OUTPUT_FILE>
```
## Command-Line Arguments

-v, --vcf <VCF_FILES>: Input VCF files. This can be multiple VCF files separated by spaces.

-l, --vcf-list <VCF_LIST_FILE>: File containing a list of VCF files, one per line. This option is mutually exclusive with -v.

-r, --reference <REFERENCE_FASTA>: Reference FASTA file (required).

-o, --output <OUTPUT_FILE>: Output file name for the distance matrix (required).

-w, --workers <NUM_WORKERS>: Number of worker threads to use (optional, defaults to the number of logical CPUs).

## Example

To calculate Fst values using multiple VCF files:
```
fstigo -v sample1.vcf sample2.vcf -r reference.fasta -o output.txt
```
Or, using a list of VCF files:
```
fstigo -l vcf_list.txt -r reference.fasta -o output.txt
```
## Output

The output is a tab-delimited file containing the pairwise Fst values for all samples provided. The diagonal of the matrix represents zero, indicating no genetic difference between a sample and itself.

## Multi-Threading

FSTiGO uses Rayon to perform parallel computations. You can specify the number of worker threads using the -w option, or it will default to the number of logical CPUs available.

## Dependencies

FSTiGO relies on the following Rust libraries:

clap: For parsing command-line arguments.

indicatif: For displaying a progress bar.

rayon: For parallel processing.

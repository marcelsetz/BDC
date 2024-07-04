#!/usr/bin/env python3

"""
This script processes ILLUMINA fastq files to calculate average quality scores per position.

It supports two modes:
1. Chunk mode: Splits the data into chunks, calculates the sum and count of the positions,
   and averages the quality scores.
2. Combine mode: Calculates the total mean quality scores of the entire file.

Usage:
    python script_name.py --chunk [--combine] [-o <output_csv_file>] <fastq_file1> <fastq_file2> ...

Arguments:
    --chunk             Run the program in chunk mode; get the sum and count of the positions.
    --combine           Run the program in combine mode; calculating the total mean of a file.
    -o                  Optional CSV file to save the output. If not provided, output is printed to stdout.
    fastq_files         At least one ILLUMINA fastq file to process.

Example:
    python script_name.py --chunk -o output.csv sample1.fastq sample2.fastq
    python script_name.py --combine -o output.csv sample1.fastq
"""

# METADATA
__author__ = "Marcel Setz"
__version__ = 2.0

# IMPORTS
import argparse
import csv
import sys
from pathlib import Path
import multiprocessing as mp


def chunks(number, mysize):
    """Returns the chunks."""
    mychunks = []
    for i in range(mysize):
        start = int(i * len(number) / mysize)
        end = int((i + 1) * len(number) / mysize)
        mychunks.append(number[start:end])

    return mychunks


def read_fastq(fastq_file):
    """Reads the files"""
    quality = []
    qual = True

    with open(fastq_file, encoding='UTF-8') as fastq:
        while qual:
            header = fastq.readline()
            nucleotides = fastq.readline()
            strand = fastq.readline()
            qual = fastq.readline().rstrip()
            if qual:
                quality.append(qual)
            if header or nucleotides or strand:
                pass

    return quality


def calculate_quals(quality):
    """Calculates quality scores"""
    results = []
    for qual in quality:
        for item, checker in enumerate(qual):
            try:
                results[item] += ord(checker) - 33
            except IndexError:
                results.append(ord(checker) - 33)
    return results


def generate_output(average_phredscores, output_file):
    """Generates the output for the file('s)"""
    if output_file is None:
        csv_writer = csv.writer(sys.stdout, delimiter=',')
    else:
        csv_writer = csv.writer(output_file.open('w', encoding='UTF-8', newline=''), delimiter=',')

    for i, score in enumerate(average_phredscores):
        csv_writer.writerow([i, score])


def main():
    # Parse arguments
    parser = argparse.ArgumentParser(
        description="Script for Assignment 3 of the Big Data Computing course."
    )

    # Create the mutually exclusive group for the mode
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--chunk",
        action="store_true",
        help="Run the program in chunk mode; get the sum and count of the positions",
    )
    mode.add_argument(
        "--combine",
        action="store_true",
        help="Run the program in combine mode; calculating the total mean of a file",
    )

    # Create group with combine mode arguments
    combine_args = parser.add_argument_group(title="Arguments when ran in combine mode")
    combine_args.add_argument(
        "-o",
        action="store",
        dest="output_file",
        type=Path,
        required=False,
        help="CSV file output should be saved to. Default is to write output to STDOUT."
    )

    parser.add_argument("fastq_files", nargs='+', help="At least 1 ILLUMINA fastq file to process")
    args = parser.parse_args()

    for file in args.fastq_files:
        qualities = read_fastq(file)
        if args.chunk:
            qual_chunked = chunks(qualities, 4)
            with mp.Pool(mp.cpu_count()) as pool:
                phredscores = pool.map(calculate_quals, qual_chunked)

            phredscores_sum = [sum(i) for i in zip(*phredscores)]
            counts = [len(chunk) for chunk in qual_chunked]
            combined_counts = sum(counts)
            phredscores_avg = [score / combined_counts for score in phredscores_sum]

            if len(args.fastq_files) > 1:
                output_file = args.output_file.with_name(f"{file}.{args.output_file.name}") if args.output_file else None
            else:
                output_file = args.output_file

            generate_output(phredscores_avg, output_file)

        elif args.combine:
            phredscores = calculate_quals(qualities)
            average_phredscores = [score / len(qualities) for score in phredscores]

            output_file = args.output_file
            generate_output(average_phredscores, output_file)


if __name__ == "__main__":
    main()

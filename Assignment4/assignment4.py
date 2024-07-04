#!/usr/bin/env python3

"""
This script processes ILLUMINA fastq files using MPI for parallel computing to calculate average quality scores per position.
It reads chunks of fastq files, distributes tasks to workers, calculates PHRED scores, and outputs the results either to stdout or an optional CSV file.

The script uses mpi4py for MPI communication and argparse for command-line argument parsing.

Usage:
    python assignment4.py -o <output_csv_file> <fastq_file1> <fastq_file2> ...

Arguments:
    -o                  Optional CSV file to save the output. Default is to write output to STDOUT.
    fastq_files         At least one Illumina FastQ Format file to process.

Example:
    python assignment4.py -o output.csv /data/datasets/rnaseq_data/fastq_file1.fastq /data/datasets/rnaseq_data/fastq_file2.fastq
"""

# METADATA
__author__ = "Marcel Setz"
__version__ = 1.0

# IMPORTS
import argparse
from pathlib import Path
from mpi4py import MPI

def parse_args():
    """Parses the CLI arguments given to the script."""
    parser = argparse.ArgumentParser(
        description="Script for Assignment 4 of the Big Data Computing course."
    )
    parser.add_argument(
        "-o",
        dest="output_file",
        type=Path,
        required=False,
        help="CSV file output should be saved to. Default is to write output to STDOUT."
    )
    parser.add_argument(
        "fastq_files",
        type=Path,
        nargs="+",
        help="At least 1 Illumina FastQ Format file to process."
    )
    return parser.parse_args()

def read_fastq_chunk(file_path, start, stop):
    """Reads a chunk of a FastQ file and returns the quality scores."""
    quality_scores = []
    with open(file_path, 'r', encoding='UTF-8') as file:
        file.seek(start)
        while file.tell() < stop:
            quality = file.readline().strip()
            if quality:
                quality_scores.append(quality)
    return quality_scores

def calculate_phred_scores(quality_scores):
    """Calculates the sum and count of PHRED scores for a chunk."""
    phred_sums = []
    counts = []
    for line in quality_scores:
        for index, char in enumerate(line):
            score = ord(char) - 33
            if len(phred_sums) <= index:
                phred_sums.append(score)
                counts.append(1)
            else:
                phred_sums[index] += score
                counts[index] += 1
    return phred_sums, counts

def process_results(results, output_file, fastq_files):
    """Processes and outputs the final results."""
    combined_sums = {}
    combined_counts = {}

    for file, (sums, counts) in results.items():
        if file not in combined_sums:
            combined_sums[file] = sums
            combined_counts[file] = counts
        else:
            for i, sum_value in enumerate(sums):
                combined_sums[file][i] += sum_value
                combined_counts[file][i] += counts[i]

    for fastq_file in fastq_files:
        averages = [combined_sums[fastq_file][i] / combined_counts[fastq_file][i] for i in range(len(combined_sums[fastq_file]))]
        if output_file:
            output_path = output_file.parent / f"{fastq_file.name}.{output_file.name}"
            with open(output_path, "w", encoding="UTF-8") as file:
                for i, score in enumerate(averages):
                    file.write(f"{i},{score}\n")
        else:
            print(f"{file.name}:")
            for i, score in enumerate(averages):
                print(f"{i},{score}")

def main():
    args = parse_args()
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nproc = comm.Get_size()

    if rank == 0:
        chunks = []
        for file_path in args.fastq_files:
            file_size = file_path.stat().st_size
            chunk_size = file_size // nproc
            for i in range(nproc):
                start = i * chunk_size
                stop = start + chunk_size if i < nproc - 1 else file_size
                chunks.append((file_path, start, stop))
    else:
        chunks = None

    chunk = comm.scatter(chunks, root=0)
    quality_scores = read_fastq_chunk(*chunk)
    phred_sums, counts = calculate_phred_scores(quality_scores)
    result = (chunk[0], (phred_sums, counts))

    all_results = comm.gather(result, root=0)

    if rank == 0:
        results = {}
        for res in all_results:
            if res[0] not in results:
                results[res[0]] = res[1]
            else:
                for i in range(len(res[1][0])):
                    results[res[0]][0][i] += res[1][0][i]
                    results[res[0]][1][i] += res[1][1][i]
        process_results(results, args.output_file, args.fastq_files)

if __name__ == "__main__":
    main()

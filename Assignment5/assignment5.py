#!/usr/bin/env python3

"""
Module for analyzing genomic features from GBFF files using PySpark.

This module contains functions to initialize a Spark session,
extract and process genomic feature data, and perform various analyses
on the processed data.
"""

# METADATA
__author__ = "Marcel Setz"
__version__ = 1.0

# IMPORTS
from pyspark.sql import Row, SparkSession
from pyspark.sql.functions import col, explode, split

# GLOBALS
CODING_KEYS = ["CDS"]
NON_CODING_KEYS = ["ncRNA", "rRNA", "gene"]
ALL_DESIRED_KEYS = [*CODING_KEYS, *NON_CODING_KEYS]

def initialize_spark_session():
    """
    Initializes and returns a Spark session with the necessary configurations.
    
    Returns:
        SparkSession: Configured Spark session.
    """
    spark = SparkSession.builder \
        .appName("GenomicFeatureAnalysis") \
        .master("local[4]") \
        .config("spark.executor.memory", "8g") \
        .config("spark.driver.memory", "8g") \
        .config("spark.jars.packages", "com.databricks:spark-xml_2.12:0.15.0") \
        .getOrCreate()
    return spark

def extract_feature_data(record):
    """
    Extracts feature information from a single record in a GBFF file.
    
    Args:
        record: A single record from the GBFF file.
        
    Returns:
        Row: A Row object containing the identifier, organism, and features.
    """
    record_str = record.value
    row = {}
    features = []
    at_features = False
    for line in record_str.split("\n"):
        if at_features:
            if not line.startswith("     "):
                break
            if line[5] != " ":
                features.append(line.strip().split())
            continue

        if line.startswith("LOCUS"):
            row["identifier"] = line.split()[1]
        elif line.startswith("SOURCE"):
            row["organism"] = line.strip().split(maxsplit=1)[1]
        elif line.startswith("FEATURES"):
            at_features = True

    row["features"] = "//".join(
        [f"{i}, {f[0]}, {f[1]}" for i, f in enumerate(features, start=1)]
    )
    return Row(**row)

def decompose_feature_info(row):
    """
    Decomposes the feature-specific info from a feature row into separate columns.
    
    Args:
        row: A Row object containing features data.
        
    Returns:
        Row: A Row object with separated feature information.
    """
    str_feat = row.features
    info = str_feat.split(", ")
    new_row = Row(
        identifier=row.identifier,
        organism=row.organism,
        feature_index=int(info[0]),
        key=info[1],
        location=info[2]
    )
    return new_row

def parse_feature_location(row):
    """
    Parses the location string of a feature into separate columns.
    
    Args:
        row: A Row object containing the feature location.
        
    Returns:
        Row: A Row object with separated location information.
    """
    location = row.location
    complement = False
    if location.startswith("complement"):
        complement = True
        location = location[11:-1]
    start, stop = location.split("..")
    new_row = Row(
        identifier=row.identifier,
        organism=row.organism,
        feature_index=row.feature_index,
        key=row.key,
        start=int(start),
        stop=int(stop),
        complement=complement,
    )
    return new_row

def create_feature_dataframe(spark, file):
    """
    Creates a DataFrame with features extracted from a .gbff file.
    
    Args:
        spark: SparkSession object.
        file: Path to the .gbff file.
        
    Returns:
        DataFrame: DataFrame containing the extracted features.
    """
    df_features = spark.read.text(file, lineSep="//\n")
    extracted_rdd = df_features.rdd.map(extract_feature_data)
    extracted_df = extracted_rdd.toDF()
    exploded_df = extracted_df.withColumn("features", explode(split("features", "//")))
    split_features_rdd = exploded_df.rdd.map(decompose_feature_info)
    split_features_df = split_features_rdd.toDF()
    filtered_df = split_features_df.filter(
        col("key").isin(ALL_DESIRED_KEYS) &
        col("location").rlike(r"^(?:complement\()?\d+\.{2}\d+\)?$")
    )
    final_rdd = filtered_df.rdd.map(parse_feature_location)
    final_df = final_rdd.toDF()
    return final_df

def exclude_coding_genes(features_df):
    """
    Excludes gene features that have a CDS feature within or spanning their location.
    
    Args:
        features_df: DataFrame containing the features.
        
    Returns:
        DataFrame excluding certain gene features.
    """
    all_genes = features_df.filter(col("key").like("gene")).alias("gene")
    all_cds = features_df.filter(col("key").like("CDS")).alias("cds")

    join_condition = [
        col("gene.identifier") == col("cds.identifier"),
        col("gene.start") <= col("cds.start"),
        col("gene.stop") >= col("cds.stop"),
        col("gene.complement") == col("cds.complement"),
    ]
    only_non_coding_genes = all_genes.join(all_cds, on=join_condition, how="left_anti")
    return features_df.filter(~col("key").like("gene")).union(only_non_coding_genes)

def question_one(df_features):
    """Calculate the average number of features per Archaea genome."""
    return df_features.groupBy("identifier").count().agg({"count": "mean"}).first()[0]


def question_two(df_features):
    """Calculate the ratio of coding to non-coding features."""
    coding_count = df_features.filter(df_features["key"].isin(CODING_KEYS)).count()
    noncoding_count = df_features.filter(df_features["key"].isin(NON_CODING_KEYS)).count()
    return coding_count / noncoding_count if noncoding_count != 0 else float('inf')


def question_three(df_features):
    """Calculate the minimum and maximum number of proteins for all organisms."""
    protein_counts = (
        df_features
        .filter(df_features.key.isin(CODING_KEYS))
        .groupBy("organism")
        .agg({"key": "count"})
        .select("count(key)")
        .sort("count(key)")
        .collect()
    )
    return protein_counts[0][0], protein_counts[-1][0]


def question_four(df_features, output_path):
    """Remove all non-coding features and save the resulting DataFrame as a Parquet file."""
    coding_df = df_features.filter(df_features["key"] == "CDS")
    coding_df.write.parquet(output_path)


def question_five(df_features):
    """Calculate the average length of a feature."""
    df_with_length = df_features.withColumn("length", col("stop") - col("start"))
    return df_with_length.agg({"length": "mean"}).first()[0]


def main():
    """
    Main function to initialize Spark, create a feature DataFrame,
    and perform various analyses on the genomic features.
    """

    spark = initialize_spark_session()
    file_path = "/data/datasets/NCBI/refseq/ftp.ncbi.nlm.nih.gov/refseq/release/archaea/archaea.1.genomic.gbff"
    output = "output/features.parquet"

    df_features = create_feature_dataframe(spark, str(file_path))
    df_features = exclude_coding_genes(df_features)

    avg_features = question_one(df_features)
    print(f"Average number of features per Archaea genome: {avg_features:.2f}")

    ratio = question_two(df_features)
    print(f"Ratio of coding to non-coding features: {ratio:.2f}")

    min_proteins, max_proteins = question_three(df_features)
    print(f"Minimum number of proteins: {min_proteins}, Maximum number of proteins: {max_proteins}")

    question_four(df_features, output)
    print(f"Coding features saved to {output}")

    avg_length = question_five(df_features)
    print(f"Average length of a feature: {avg_length:.2f}")

if __name__ == "__main__":
    main()

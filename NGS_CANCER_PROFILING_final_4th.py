# Databricks notebook source
# MAGIC %md #Setting Up NGS Platform and Configuring

# COMMAND ----------

pip install BIo

# COMMAND ----------

import pyspark
from Bio import SeqIO
from pyspark.sql import SparkSession
from pyspark.sql import Row

# COMMAND ----------

spark.conf.set(f"fs.azure.sas.raw.ngddatacancerprofiling.blob.core.windows.net",'?sv=2022-11-02&ss=bfqt&srt=sco&sp=rwdlacupiytfx&se=2023-08-04T01:56:38Z&st=2023-08-03T17:56:38Z&spr=https&sig=VqJDiNMv2jVQiKpLtkI82Bgf4cAeKAAfZ943kI5GyOQ%3D')

# COMMAND ----------

dbutils.fs.ls('wasbs://raw@ngddatacancerprofiling.blob.core.windows.net')

# COMMAND ----------

# dbutils.fs.mount(
#   source='wasbs://raw@ngddatacancerprofiling.blob.core.windows.net',
#   mount_point='/mnt/raw',
#   extra_configs={'fs.azure.account.key.ngddatacancerprofiling.blob.core.windows.net':'Wu6/hgfQirWHysHsOwD3mbj7Jl0ewxFyd0E6pfKD5gqZfbhQictzDGRw4IzDY/crFiuh6wRcBAV8+AStmN7i9Q=='}
# )

# COMMAND ----------

# MAGIC %md #Reading cancer Fastq files and Reference genome

# COMMAND ----------

def read_fastq_as_rdd_and_dataframe(file_path):
    # Create a SparkSession
    spark = SparkSession.builder \
        .appName("FASTQ Reader") \
        .getOrCreate()

    # Read FASTQ file using Biopython
    fastq_records = SeqIO.parse(file_path, "fastq")
    # for record in fastq_records:
    #     print("ID:", record.id)
    #     print("Sequence:", record.seq)
    #     print("Quality:", record.letter_annotations["phred_quality"])
    #     print("---")
    record_dicts = []
    for record in fastq_records:
        record_dict = {
            "id": record.id,
            "sequence": str(record.seq),
            "quality": str(record.letter_annotations["phred_quality"])
        }
        record_dicts.append(record_dict)
    # print(record_dicts)
    # Create an RDD from the list of dictionaries
    rdd = spark.sparkContext.parallelize(record_dicts)
    # rdd.take(5).foreach(print)
    # Convert RDD to DataFrame
    # df = rdd.map(lambda x: Row(**x)).toDF()
    return rdd 

# COMMAND ----------

# file_path = "/dbfs/mnt/raw/test_ngs.fastq"
# # file_path='wasbs://raw@ngddatacancerprofiling.blob.core.windows.net/test_ngs.fastq'
# cancer_rdd=read_fastq_as_rdd_and_dataframe(file_path)

# COMMAND ----------

# for row in cancer_rdd.collect():
#     print(row)

# COMMAND ----------

# df=cancer_rdd.map(lambda x: Row(**x)).toDF()

# COMMAND ----------

# df.show()

# COMMAND ----------

file_path = "/dbfs/mnt/raw/partaa.fastq"
# file_path='wasbs://raw@ngddatacancerprofiling.blob.core.windows.net/test_ngs.fastq'
cancer_rdd=read_fastq_as_rdd_and_dataframe(file_path)
df_cancer=cancer_rdd.map(lambda x: Row(**x)).toDF()

# COMMAND ----------

# for row in cancer_rdd.collect():
#     print(row)

# COMMAND ----------


# df_cancer.show()

# COMMAND ----------

# df.count()

# COMMAND ----------

from pyspark import SparkContext, SparkConf
from pyspark.sql import SparkSession

# Configure Spark
spark = SparkSession.builder \
        .appName("FASTQ Reader") \
        .getOrCreate()

# Replace 'file:///path/to/your/GCF_000001405.40_GRCh38.p14_genomic.fna' with the actual file path
genome_file_path = '/mnt/raw/GCF_000001405.40_GRCh38.p14_genomic.fna'

# Load the genome file into an RDD
genome_rdd = spark.read.text(genome_file_path)

# COMMAND ----------


# genome_rdd,df

# bam_df = fastq_df.join(reference_genome_df, on="sequence")

# bam_df = df.join(genome_rdd,df, on="sequence")

# COMMAND ----------

# genome_rdd.show()

# COMMAND ----------

# genome_rdd.count()

# COMMAND ----------

# MAGIC %md #Cleaning of Dataframe created from reference genome

# COMMAND ----------

# Create a Spark session
spark = SparkSession.builder.appName("Cleaning Reference Genome").getOrCreate()

# Load the reference genome
reference_genome_df = spark.read.format("fasta").load("data/reference_genome.fasta")

# Remove repetitive sequences
repetitive_sequences = reference_genome_df.filter(reference_genome_df["sequence"].rlike("N+"))
reference_genome_df = reference_genome_df.subtract(repetitive_sequences)

# Remove low-quality sequences
low_quality_sequences = reference_genome_df.filter(reference_genome_df["quality"].lt(20))
reference_genome_df = reference_genome_df.subtract(low_quality_sequences)

# Trim adapters
adapters = reference_genome_df.filter(reference_genome_df["sequence"].startswith("AGATCGGAAGAGCACACGTC"))
reference_genome_df = reference_genome_df.subtract(adapters)

# Fill gaps
gaps = reference_genome_df.filter(reference_genome_df["sequence"].contains("-"))
reference_genome_df = reference_genome_df.subtract(gaps)

# Save the cleaned reference genome
reference_genome_df.write.save("data/cleaned_reference_genome.fasta", format="fasta")

# COMMAND ----------

# genome_rdd['value']=genome_rdd['value'].apply(str.upper)
# import pyspark.sql.functions as F

# df_reference=genome_rdd.select(F.upper(genome_rdd.value))


# COMMAND ----------

genome_path = '/mnt/raw/GCF_000001405.40_GRCh38.p14_genomic.fna'
fastq_path = '/mnt/raw/ERR10360601_1.fastq'

# Load FASTQ data as DataFrame (assuming it's in FASTQ format)
fastq_df = spark.read.format('fastq').load(fastq_path)

# Load reference genome as RDD (assuming it's in FASTA format)
genome_rdd = spark.sparkContext.textFile(genome_path)

# Step 3: Preprocess the data (optional)
# Perform preprocessing steps as required, e.g., quality filtering, adapter trimming, etc.
# You can use PySpark functions or external tools for this step.

# Step 4: Perform read alignment using an external aligner (e.g., Bowtie, BWA, or HISAT2)
# Use Python's subprocess module to call the aligner from the shell
import subprocess

# Replace 'path/to/aligner' with the actual path to the aligner binary on Databricks
# Replace 'path/to/output' with the desired output path for alignment results (SAM/BAM file)
aligner_command = ['/dbfs/mnt/raw', '--options', genome_path, fastq_path, '-S', '/dbfs/mnt/raw']
subprocess.run(aligner_command, shell=True)

# Step 5: Process and analyze the alignment results
# Load the alignment results (SAM/BAM file) into a DataFrame (assuming it's in SAM format)
alignment_results_path = '/dbfs/path/to/output'
alignment_df = spark.read.format('sam').load(alignment_results_path)

# Perform any necessary processing and analysis on the alignment DataFrame
# For example, count mapped reads, calculate alignment quality scores, etc.

# Display or save the results as required
alignment_df.show()

# COMMAND ----------

df_reference.show()


# COMMAND ----------

df_cancer.show()

# COMMAND ----------

df_reference.count()

# COMMAND ----------

df_cancer.count()

# COMMAND ----------

repetitive_sequences = df_reference.filter(df_reference["upper(value)"].rlike("N"))
reference_genome_df = df_reference.subtract(repetitive_sequences)

# COMMAND ----------

reference_genome_df.count()

# COMMAND ----------

reference_genome_df.show()

# COMMAND ----------

df_cancer.show()

# COMMAND ----------

pip install subprocess.run

# COMMAND ----------

# MAGIC %md #  Alignment

# COMMAND ----------

import subprocess

# COMMAND ----------

genome_path = '/mnt/raw/GCF_000001405.40_GRCh38.p14_genomic.fna'
fastq_path = '/mnt/raw/partaa.fastq'

# COMMAND ----------

aligner_command = ['/dbfs/mnt/raw/bowtie-1.3.1-macos-x86_64/bowtie', '--options', genome_path, fastq_path, '-S', '/dbfs/mnt/raw/output']
subprocess.run(aligner_command, shell=True)

# COMMAND ----------

# MAGIC %md #Variant Calling

# COMMAND ----------

pip install hail

# COMMAND ----------

import hail as hl
hl.init()

# COMMAND ----------

input_bam_path = '/dbfs/path/to/aligned_reads.bam'
mt = hl.import_bam(input_bam_path, reference_genome='GRCh38')

# COMMAND ----------

mt_snps = mt.filter_rows(mt.alleles.is_snp())

# Show SNPs (You can perform additional analyses or write the results to a file)
mt_snps.rows().show()

# COMMAND ----------

# Filter the matrix table to keep only insertions
mt_insertions = mt.filter_rows(mt.alleles.is_insertion())

# Show insertions (You can perform additional analyses or write the results to a file)
mt_insertions.rows().show()

# COMMAND ----------

mt_deletions = mt.filter_rows(mt.alleles.is_deletion())

# Show deletions (You can perform additional analyses or write the results to a file)
mt_deletions.rows().show()





# COMMAND ----------

# Call structural variations using hl.lumpy
lumpy_result = hl.lumpy(mt)

# Show structural variations (You can perform additional analyses or write the results to a file)
lumpy_result.show()

# COMMAND ----------

# MAGIC %md #Gene expression

# COMMAND ----------

import pandas as pd
import numpy as np
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri

# Install the required R packages (DESeq2 and edgeR) using rpy2
r_install_packages = """
install.packages("DESeq2")
install.packages("edgeR")
"""

robjects.r(r_install_packages)

def differential_gene_expression_analysis(data, condition_col, method='DESeq2'):
    """
    Perform differential gene expression analysis using DESeq2 or edgeR.

    Parameters:
        data (pandas.DataFrame): DataFrame containing gene expression data. Rows represent genes, columns represent samples.
        condition_col (str): The name of the column in the DataFrame that represents the experimental condition or group for each sample.
        method (str): The method to use for differential gene expression analysis. Options: 'DESeq2' (default) or 'edgeR'.

    Returns:
        pandas.DataFrame: DataFrame containing the differential gene expression results.
    """
    # Ensure that rpy2 will convert pandas DataFrame to R DataFrame
    pandas2ri.activate()

    # Load the required R packages (DESeq2 and edgeR)
    robjects.r(f"library({method})")

    # Convert the pandas DataFrame to an R DataFrame
    r_data = pandas2ri.py2ri(data)

    # Perform normalization and differential expression analysis
    if method == 'DESeq2':
        # Using DESeq2 for differential expression analysis
        robjects.r(f"dds <- DESeqDataSetFromMatrix(countData=r_data, colData=as.data.frame(r_data['{condition_col}']), design=~{condition_col})")
        robjects.r("dds <- DESeq(dds)")
        result = robjects.r("res <- results(dds)")

    elif method == 'edgeR':
        # Using edgeR for differential expression analysis
        robjects.r(f"y <- DGEList(counts=r_data, group=as.factor(r_data['{condition_col}']))")
        robjects.r("y <- calcNormFactors(y)")
        robjects.r("design <- model.matrix(~as.factor(r_data['{condition_col}']))")
        robjects.r("y <- estimateDisp(y, design)")
        robjects.r("fit <- glmFit(y, design)")
        result = robjects.r("res <- glmLRT(fit)")

    # Convert the R DataFrame result back to a pandas DataFrame
    result_df = pandas2ri.ri2py_dataframe(result)

    return result_df
Example usage:


# Assuming you have loaded your gene expression data into a pandas DataFrame 'gene_expression_df'

# Perform differential gene expression analysis using DESeq2
deseq2_results = differential_gene_expression_analysis(gene_expression_df, condition_col='Condition', method='DESeq2')

# Perform differential gene expression analysis using edgeR
edger_results = differential_gene_expression_analysis(gene_expression_df, condition_col='Condition', method='edgeR')


# COMMAND ----------

#Perform differential gene expression analysis using DESeq2
deseq2_results = differential_gene_expression_analysis(gene_expression_df, condition_col='Condition', method='DESeq2')

# Perform differential gene expression analysis using edgeR
edger_results = differential_gene_expression_analysis(gene_expression_df, condition_col='Condition', method='edgeR')

import gzip
import os
import re
import subprocess
import joblib
import pandas as pd 
import dask.dataframe as dd
import numpy as np
from pathlib import Path

import scipy
import scipy.sparse as sp

from sklearn.cluster import KMeans
from sklearn.preprocessing import Normalizer
from sklearn.feature_extraction.text import HashingVectorizer
from sklearn.feature_extraction.text import TfidfVectorizer
directory_path = Path(__file__).parent.absolute()

VCF_BASE_URL = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/"
VCF_DIR = str(directory_path / "vcf_files")
TSV_DIR = str(directory_path / "tsv_files")

CLUSTER_FILE = str(directory_path / "cluster.pkl")
VECTOR_FILE = str(directory_path / "vector.npz")
LABELS_FILE = str(directory_path / "labels.npy")

CHUNK_SIZE = 65536

base_columns = ["CHROM", "POS", "REF", "ALT", "GENOME"] 
legend_columns = ["url", "md5", "Data collection", "Data type", "Analysis group", "Sample", "Population", "Data reuse policy"]
number_chromesomes = 23
chromosomes = range(1, number_chromesomes)

vectorizer = HashingVectorizer(analyzer="char", ngram_range=(3, 3), n_features=9)

file_suffix = ".tsv"

pattern = r"([01]\|[01](?:\t[01]\|[01])*)"

TESTING = True

if TESTING:
    chromosomes = [14, 20]

def create_data_storage_folders():
    os.makedirs(VCF_DIR, exist_ok=True)

    # os.makedirs(TSV_DIR, exist_ok=True)

def get_VCF_uri(chromosome):
    #        ALL.chr14          .shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
    return f"{VCF_BASE_URL}ALL.chr{chromosome}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"

def get_VCF_path(chromosome):

    return f"{VCF_DIR}{os.sep}Chr{chromosome}.vcf.gz"

def download_VCF_file(chromosome):
    
    VCF_PATH = get_VCF_path(chromosome)
    DOWNLOAD_URI = get_VCF_uri(chromosome)

    cmd = f'curl --ftp-pasv -o "{VCF_PATH}"  "{DOWNLOAD_URI}"'

    subprocess.run(cmd, shell=True, check=True)

def delete_VCF_file(chromosome):
    file_path = get_VCF_path(chromosome)
    os.remove(file_path)

def process_vcf_in_chunks(file_path, pattern, individuals, chunk_size=CHUNK_SIZE, npartitions=8):
    df_main = None
    i = 0
    
    with gzip.open(file_path, 'rt') as f:
        chunk = []
        for line in f:
            if line.startswith("#"):
                continue  # Ignore comments
            
            chunk.append(line.strip())

            if sum(len(l) for l in chunk) > chunk_size:

                matches = re.findall(pattern, "\n".join(chunk))
                
                if matches:
                    df_chunk = dd.from_pandas(pd.DataFrame(matches, columns=["Genotype"]), npartitions=npartitions)
                    
                    if df_main is None:
                        df_main = df_chunk
                    else:
                        df_main = dd.concat([df_main, df_chunk])

                chunk = [] 
                i += 1

        if chunk:
            matches = re.findall(pattern, "\n".join(chunk))
            
            if matches:
                df_chunk = dd.from_pandas(pd.DataFrame(matches, columns=["Genotype"]), npartitions=npartitions)
                
                if df_main is None:
                    df_main = df_chunk
                else:
                    df_chunk = df_chunk.compute()
                    df_main.append(df_chunk)

    df_main.columns = individuals

    return df_main


def load_data(file, individuals):
    ddf = dd.read_csv(file, sep="\t", header=None, dtype=str, blocksize=CHUNK_SIZE)

    #The fourth column no longer is needed -> was originally genome, but that gets split off into individuals
    columns = base_columns[0:4] + individuals
    ddf.columns = columns

    return ddf

def pad_missing_individuals(df, expected_individuals):
    
    current_individuals = df.columns

    missing_individuals = set(expected_individuals) - set(current_individuals)
    for individual in missing_individuals:
        df[individual] = '0'  
    
    df = df[expected_individuals]

    return df

def remove_pipe_characters(df):
    df = df.fillna('') 
    return df.map(lambda x: x.replace('|', '') if isinstance(x, str) else x)


def stream_process_vcf_file(file_path, pattern, individual_vectors, individuals, vectorizer, chunk_size=CHUNK_SIZE):
    tail_buffers = {ind: "" for ind in individuals}
    
    with gzip.open(file_path, 'rt') as f:
        chunk = []
        for line in f:
            if line.startswith("#"):
                continue

            chunk.append(line.strip())

            if sum(len(l) for l in chunk) > chunk_size:
                tail_buffers = process_chunk(chunk, pattern, individual_vectors, individuals, vectorizer, tail_buffers)
                chunk = []

        if chunk:
            tail_buffers = process_chunk(chunk, pattern, individual_vectors, individuals, vectorizer, tail_buffers)

def process_chunk(chunk, pattern, individual_vectors, individuals, vectorizer, tail_buffers):
    matches = re.findall(pattern, "\n".join(chunk))
    if not matches:
        return tail_buffers

    for match in matches:
        genotypes = match.split("\t")
        for idx, individual in enumerate(individuals):
            if individual not in individual_vectors:
                individual_vectors[individual] = sp.csr_matrix((1, vectorizer.n_features))

            # Remove '|' and build full sequence with tail from last chunk
            sequence = tail_buffers[individual] + genotypes[idx].replace("|", "")
            tail_buffers[individual] = sequence[-2:]  # retain 2 for 3-mers

            kmers = " ".join(sequence[i:i+3] for i in range(len(sequence) - 2))
            if kmers:
                vec = vectorizer.transform([kmers])
                individual_vectors[individual] += vec

    return tail_buffers


def process_dask_dataframe(ddf, individual_vectors):
    ddf = ddf.map_partitions(remove_pipe_characters)

    vectorizer = HashingVectorizer(analyzer="char", ngram_range=(3, 3), n_features=9)  

    for individual in ddf.columns:
        if individual not in individual_vectors:
            individual_vectors[individual] = sp.csr_matrix((1, 9)) 

        batches = ddf[individual].to_delayed()
        for delayed_batch in batches:
            batch = delayed_batch.compute()

            kmer_sequence = " ".join(batch.astype(str))
            
            vectorized = vectorizer.transform([kmer_sequence])  

            individual_vectors[individual] += vectorized

    return individual_vectors


def train_clustering_model(n_clusters=5, save_results=True):
    global individual_vectors       
    
    if individual_vectors is None:
        raise ValueError("No data found in individual_vectors. Run update_hash_vectorizer first!")

    normalizer = Normalizer()
    individual_vectors_norm = normalizer.fit_transform(individual_vectors)

    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    clusters = kmeans.fit_predict(individual_vectors_norm)

    if save_results:
        save_clustering_results(kmeans, clusters)
        print("Saved clustering results")

    return clusters, kmeans

def save_sparse_matrix(filename, matrix):
    scipy.sparse.save_npz(filename, matrix)

def load_sparse_matrix(filename):
    return scipy.sparse.load_npz(filename)

def save_clustering_results(model, clusters):
    joblib.dump(model, CLUSTER_FILE)
    np.save(LABELS_FILE, clusters)

def load_clustering_results():
    model = joblib.load(CLUSTER_FILE)
    clusters = np.load(LABELS_FILE)
    return model, clusters

create_data_storage_folders()

legend_df = pd.read_csv(str(directory_path) + f"{os.sep}Data{os.sep}Key" + file_suffix, delimiter="\t")
legend_df.columns = legend_columns

sample_entries = legend_df['Sample'].array

individuals = []
for v in sample_entries:
    split = re.split(r",\s?", v)
    for sample in split:
        if not sample in individuals:
            individuals.append(sample)

individual_vectors = {}
individual_index = {}

for chromosome in chromosomes:
    print(f"Beginning Chromosome {chromosome}")

    genome_column_labels = legend_df.loc[chromosome - 1, "Sample"].split(",")

    if not TESTING:
        print(f"Downloading chromosome {chromosome} data")
        download_VCF_file(chromosome)


    print(f"Loading chromosome {chromosome} data into dataframe")
    # ddf = load_data(get_TSV_path(chromosome), genome_column_labels)
    # ddf = ddf.iloc[:, 4:] 
    vcf_file = get_VCF_path(chromosome)

    # ddf = process_vcf_in_chunks(vcf_file, pattern, genome_column_labels)

    # ddf = pad_missing_individuals(ddf, individuals)

    stream_process_vcf_file(vcf_file, pattern, individual_vectors, individuals, vectorizer)


    # individual_vectors = process_dask_dataframe(ddf, individual_vectors)

    print("Deleting files")

    if not TESTING:
        delete_VCF_file(chromosome)
    

individual_vectors = sp.vstack([individual_vectors[ind] for ind in sorted(individual_vectors.keys())])

save_sparse_matrix(VECTOR_FILE, individual_vectors)
print("Saved the HashVectorizer results")

clusters, model = train_clustering_model(n_clusters=2)

for i, cluster in enumerate(clusters):
    print(f"Individual {i} -> Cluster {cluster}")
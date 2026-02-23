import pandas as pd
import os

def split_csv(file_path, chunk_size=1000):
    """Splits a CSV file into smaller chunks.

    Args:
        file_path (str): Path to the CSV file.
        chunk_size (int): Number of rows per chunk.
    """
    file_name, file_extension = os.path.splitext(os.path.basename(file_path))
    
    reader = pd.read_csv(file_path, chunksize=chunk_size)
    
    chunk_number = 1
    for chunk in reader:
        chunk.to_csv(f"{file_name}_chunk_{chunk_number}{file_extension}", index=False)
        chunk_number += 1
        
# Example usage:
file_path = "/nesi/nobackup/uc04105/new_databases_May/GTDB_226/Bacteria_GTDB226_protein_May92025.tsv"
split_csv(file_path, chunk_size=7070931)

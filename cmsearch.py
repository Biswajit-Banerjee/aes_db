import subprocess
from tqdm.auto import tqdm
import os
import re
import pandas as pd
from multiprocessing import Pool, cpu_count

FASTA_PATH = "/home/sumon/workspace/git_repos/aes_db/data/bac_ribosome.fasta"

def cm_search(pair, fasta_path=FASTA_PATH):
    (cm_path, out_file) = pair
    command = ["/home/R2DT1_4/infernal/infernal-1.1.5/src/cmsearch", "--tblout", out_file, cm_path, fasta_path]
    result = subprocess.run(command, 
                            capture_output=True, text=True)
    with open(out_file.replace(".tblout", ".out"), 'w') as f:
        f.writelines(result.stdout)
    return out_file


if __name__ == "__main__":
    base_cm_path = "/home/sumon/workspace/git_repos/aes_db/data/all_cms"
    base_blout = "/home/sumon/workspace/git_repos/aes_db/data/all_tsv"
    os.makedirs(base_blout, exist_ok=True)
    
    # Create a list of all pairs to compare
    all_pairs = [
        (os.path.join(base_cm_path, cm), os.path.join(base_blout, cm+".tblout")) 
        for cm in os.listdir(base_cm_path) if cm.endswith(".cm")
    ]
    
    # Use multiprocessing to compare models in parallel
    with Pool(processes=cpu_count()) as pool:
        results = list(tqdm(pool.imap(cm_search, all_pairs), total=len(all_pairs)))
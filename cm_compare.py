import subprocess
from tqdm.auto import tqdm
import os
import re
import pandas as pd
from multiprocessing import Pool, cpu_count

cm_dir = "/home/sumon/repos/aes_db/data/all_cms"
all_cm_files = os.listdir(cm_dir)

def compare_models(args):
    name1, name2 = args
    cm1 = os.path.join(cm_dir, name1)
    cm2 = os.path.join(cm_dir, name2)
    result = subprocess.run(["/home/sumon/repos/aes_db/executables/cmcompare", cm1, cm2], 
                            capture_output=True, text=True)
    std_out = re.sub(" +", " ", result.stdout).split(" ")
    return (name1.replace(".cm", ""), name2.replace(".cm", ""), std_out[2])

if __name__ == "__main__":
    comparision_dt = {}

    # Create a list of all pairs to compare
    all_pairs = [(name1, name2) for name1 in all_cm_files for name2 in all_cm_files]

    # Use multiprocessing to compare models in parallel
    with Pool(processes=cpu_count()) as pool:
        results = list(tqdm(pool.imap(compare_models, all_pairs), total=len(all_pairs)))

    # Organize results into a dictionary
    for name1, name2, value in results:
        if name1 not in comparision_dt:
            comparision_dt[name1] = {}
        comparision_dt[name1][name2] = value

    # Convert to DataFrame and save as CSV
    df = pd.DataFrame(comparision_dt)
    df.to_csv("./comparision.csv")

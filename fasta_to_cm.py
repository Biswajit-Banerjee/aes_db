from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import subprocess 
import json
import os
import argparse
from datetime import datetime as dt
import multiprocessing
from functools import partial
from utils import log, modify_extension
from msa_breaker import break_msa, process_aes_df, write_stockholm
import requests

def fasta_to_stockholm(fasta_path, sto_path="", annchor_name=None):
    sto_path = modify_extension(fasta_path, "stockholm", dir_path=sto_path)
    
    records = []
    anchor_gapped_sequence = ""
    
    
    for record in SeqIO.parse(fasta_path, "fasta"):
        new_id = record.description.replace(" ", "__")
        new_record = SeqRecord(record.seq, id=new_id, description="")
        records.append(new_record)
        if annchor_name and annchor_name in record.description:
            anchor_gapped_sequence = record.seq
    
    count = SeqIO.write(records, sto_path, "stockholm")
    return sto_path, anchor_gapped_sequence

def run_cmbuild(sto_path, cm_path=""):
    cm_path = modify_extension(sto_path, "cm", dir_path=cm_path)    
    result = subprocess.run([f"cmbuild", "-F", cm_path, sto_path], 
                            capture_output=True, text=True)
    
    if result.returncode != 0:
        log(f"Error running cmbuild for {sto_path}, see below for more details", "error")
        print(result.stderr)
        return None
    
    return cm_path


def convert_json_to_dotbracket(bp_json, index_mapping=None, max_length=None):
        
    annotations = bp_json["annotations"]
    
    open_brackets = []
    closed_brackets = []
    
    for value in annotations:                
        if value["bp"] == "cWW" and value["crossing"] == "0":
            open_brackets.append(int(value['seq_id1']))
            closed_brackets.append(int(value['seq_id2']))
        
    dot_bracket = (max_length) * ["."]
    
    for symbol, index_list in [("(", open_brackets), (")", closed_brackets)]:
        for i in index_list:
            
            if index_mapping:
                mapped_i = index_mapping[i]
            else:
                mapped_i = i
            
            dot_bracket[mapped_i] = symbol
    
    dot_bracket = "".join(dot_bracket)
    
    return dot_bracket

def create_mapping(gapped_sequence, start_index=1):
    '''
    Maps anchor index to MSA index
    '''
    mapping = {}
    
    current_index = start_index
    for i, symbol in enumerate(gapped_sequence, start_index):
        if symbol != "-":
           mapping[current_index] = i
           current_index += 1
    
    with open("./mapping.json", "w") as fp:
        json.dump(mapping, fp, indent=4)
        
    return mapping

def insert_line_to_file(file_path, line, index):
    with open(file_path, "r") as f:
        contents = f.readlines()
        
    contents.insert(index, line + "\n")
    
    with open(file_path, "w") as f:
        contents = "".join(contents)
        f.write(contents)


def extract_anchor_from_msa(fasta_path, annchor_name):
    anchor_gapped_sequence = ""
    
    for record in SeqIO.parse(fasta_path, "fasta"):
        if annchor_name and annchor_name in record.description:
            anchor_gapped_sequence = record.seq
            break 
    return anchor_gapped_sequence
    

def process_fasta(fasta_file, args):
    start_time = dt.now()
    fasta_path = os.path.join(args.fasta_path, fasta_file)
    stockholm = args.stockholm
    organism  = args.organism
    pdb_id = args.pdb_id
    cm = args.cm
    chain = args.chain
    aes_csv_path = args.aes_csv_path
    
    # sto_path, gapped_anchor_seq = fasta_to_stockholm(fasta_path, stockholm, organism)
    gapped_anchor_seq = extract_anchor_from_msa(fasta_path, organism)
    log(f"Organism sequence extracted from MSA: {organism}")
    
    mapping = create_mapping(gapped_anchor_seq)
    log(f"Anchor sequence to MSA mapping created")
    
    bp_json = run_fr3d(pdb_id, chain)
    log(f"Base pairs extracted for anchor structure: {pdb_id} chain {chain}")
    
    dot_bracket = convert_json_to_dotbracket(bp_json, mapping, max_length=len(gapped_anchor_seq))
    log(f"Dot bracket created and mapped to MSA")
    
    aes_mapping = process_aes_df(aes_csv_path, mapping)
    log(f"AES mapping created")
    
    aes_records = break_msa(fasta_path, aes_mapping, dot_bracket)
    log(f"Records created for individual AES: {len(aes_records)}")
    
    all_files = write_stockholm(aes_records)
    log(f"Dot bracket added to stockholm")

    log(f"Constructing covariance model for {len(all_files)} AESs, please wait...")
    for sto_path in all_files:
        cm_path = run_cmbuild(sto_path, cm)
    
        if cm_path:
            log(f"CM file created: {cm_path}")
        else:
            log(f"Failed to create CM file for {sto_path}", "error")
            break
    
    log(f"Time taken for {fasta_file}: {dt.now() - start_time}", "debug")

    if args.clean and not args.stockholm:
        os.remove(sto_path)
        log(f"Cleaned up stockholm file: {sto_path}", "debug")
    
    return dt.now() - start_time

def run_fr3d(pdb_id, chain):
    # path, file = os.path.split(file_path)
    # output_dir = os.environ.get("OUTPUT_DIR", "./output")
    # os.makedirs(output_dir, exist_ok=True)
    # command = f"python {os.path.join(fred_path, 'classifiers/NA_pairwise_interactions.py')} --input {path} {file} -o {output_dir} -f ebi_json --chain {chain}"
    # output = subprocess.run(command.split(" "), capture_output=True)
    # log(f"FR3D: {output}", "debug")
    
    # out_file = os.path.join(output_dir, f"{file.split('.')[0]}_{chain}_basepair.json")
    url = f"https://rna.bgsu.edu/rna3dhub/rest/getChainSequenceBasePairs?pdb_id={pdb_id}&chain={chain}&only_nested=True"
    response = requests.get(url)
    return response.json()

def main():
    parser = argparse.ArgumentParser(description="Convert FASTA to Stockholm and build CM in parallel")
    parser.add_argument("fasta_path", help="Path to the dir containing input FASTA files")
    parser.add_argument("-s", "--stockholm", help="Path to store intermediate output Stockholm files (optional)")
    parser.add_argument("-c", "--cm", help="Path to store output CM files (optional)")
    parser.add_argument("--clean", action="store_true", default=False, help="Clean up intermediate Stockholm files")
    parser.add_argument("-j", "--jobs", type=int, default=multiprocessing.cpu_count(), help="Number of parallel jobs (default: number of CPU cores)")
    
    args = parser.parse_args()
    
    all_fasta_files = [file for file in os.listdir(args.fasta_path) if file.endswith((".fas", ".fasta"))]
    
    start = dt.now()
    
    with multiprocessing.Pool(processes=args.jobs) as pool:
        process_func = partial(process_fasta, args=args)
        times = pool.map(process_func, all_fasta_files)
    
    total_time = dt.now() - start
    log(f"Total time elapsed: {total_time}", "debug")

if __name__ == "__main__":
    pass
    # main()
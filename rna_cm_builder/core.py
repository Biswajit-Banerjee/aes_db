import os
import json
import requests
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio import SeqIO
from .utils import log, modify_extension
from .msa_breaker import process_aes_df, break_msa, write_stockholm, write_fasta
from .cm_builder import run_cmbuild, convert_json_to_dotbracket, create_mapping, extract_base_pairs, colocate_basepairs_in_aes, run_cm_convert
from .config import Config
from .config import Binaries as bins
from .run_ipknot import run_ipknot
from .show_ss import show_stockholm
from tqdm.auto import tqdm

def extract_anchor_from_msa(fasta_path, anchor_name):
    for record in SeqIO.parse(fasta_path, "fasta"):
        if anchor_name and anchor_name in record.description:
            print()
            return str(record.seq)
    return ""

def run_fr3d(pdb_id, chain):
    url = f"https://rna.bgsu.edu/rna3dhub/rest/getChainSequenceBasePairs?pdb_id={pdb_id}&chain={chain}&only_nested=True"
    response = requests.get(url)
    return response.json()

def process_single_msa(msa_data):
    fasta_path = msa_data['fasta_path']
    organism = msa_data['organism']
    pdb_id = msa_data['pdb_id']
    chain = msa_data['chain']
    aes_csv_path = msa_data['aes_csv_path']
    alignment_name = os.path.basename(fasta_path).split(".")[0]

    gapped_anchor_seq = extract_anchor_from_msa(fasta_path, organism)
    log(f"Organism sequence extracted from MSA: {organism}")

    mapping = create_mapping(gapped_anchor_seq)
    log(f"Anchor sequence to MSA mapping created")

    bp_json = run_fr3d(pdb_id, chain)
    log(f"Base pairs extracted for anchor structure: {pdb_id} chain {chain}")

    base_pairs = extract_base_pairs(bp_json, mapping)
    log(f"Base pairs filtered, count: {len(base_pairs)}")

    aes_mapping = process_aes_df(aes_csv_path, mapping, augment=True)
    log(f"AES mapping created, count: {len(aes_mapping)}")
    
    aes_bp_mapping = colocate_basepairs_in_aes(base_pairs, aes_mapping)
    log(f"Base pairs broken into AES chunks, count: {len(aes_bp_mapping)}")

    aes_records = break_msa(fasta_path, aes_mapping, aes_bp_mapping)
    log(f"Records created for individual AES: {len(aes_records)}")

    all_files = write_stockholm(aes_records, Config.DATA_DIR / alignment_name / "stockholms")
    log(f"Dot bracket added to stockholm")

    log(f"Constructing covariance model for {len(all_files)} AESs, please wait...")
    for sto_path in all_files: 
        cm_path  = run_cmbuild(sto_path)
        hmm_path = run_cm_convert(cm_path)
        cm1_path = run_cm_convert(cm_path, flag='-1')
        
        if cm_path and hmm_path:
            log(f"CM & HMM created: {cm_path}")
        elif cm_path:
            log(f"Failed to create HMM file for {sto_path}", "error")
        elif hmm_path:
            log(f"Failed to create CM file for {sto_path}", "error")
        else:
            log(f"Failed to create CM & HMM file for {sto_path}", "error")
            
    


def process_single_msa_using_ipknot(msa_data):
    fasta_path = msa_data['fasta_path']
    organism = msa_data['organism']
    pdb_id = msa_data['pdb_id']
    chain = msa_data['chain']
    aes_csv_path = msa_data['aes_csv_path']
    result_dir = msa_data.get("local_storage", "./results")

    gapped_anchor_seq = extract_anchor_from_msa(fasta_path, organism)
    
    mapping = create_mapping(gapped_anchor_seq)
    log(f"mapping created")
    
    bp_json = run_fr3d(pdb_id, chain)
    log(f"Base pairs extracted for anchor structure: {pdb_id} chain {chain}")

    dot_bracket = convert_json_to_dotbracket(bp_json, mapping, max_length=len(gapped_anchor_seq))
    log(f"Dot bracket created and mapped to MSA")

    aes_mapping = process_aes_df(aes_csv_path, mapping)
    log(f"AES mapping created")

    aes_records = break_msa(fasta_path, aes_mapping, dot_bracket)
    log(f"Records created for individual AES: {len(aes_records)}")

    lines = []
    for aes, [records, partial_dot_bracket] in tqdm(aes_records.items()):
        fasta_path = os.path.join(result_dir, f"aes_{aes}", "msa.fasta")
        os.makedirs(os.path.dirname(fasta_path), exist_ok=True)
        SeqIO.write(records, fasta_path, "fasta")
        sequence, conc_dot_bracket = run_ipknot(fasta_path)
        lines.append({
            "FR3D": partial_dot_bracket,
            "IPKNOT": conc_dot_bracket,
            "SEQUENCE": sequence,
            "TITLE": f"AES {aes}"
        })
        
    with open("test_1.json", "w") as f:
        json.dump(lines, f, indent=4)
            
def process_rna_data(input_file):
    with open(input_file, 'r') as f:
        data = json.load(f)

    if isinstance(data, list):  # Multiple MSA files
        with ProcessPoolExecutor() as executor:
            futures = [executor.submit(process_single_msa, msa_data) for msa_data in data]
            for future in as_completed(futures):
                future.result()
    else:  # Single MSA file
        process_single_msa(data)

    log("All processing completed successfully.")
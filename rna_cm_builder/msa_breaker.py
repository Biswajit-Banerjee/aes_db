from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
from pathlib import Path
import os
from .utils import log, modify_extension

def process_range(range_str, mapping):
    processed_ranges = []
    range_pairs = range_str.split(";")
    for pair in range_pairs:
        if not pair.strip():
            continue
        start, end = pair.split("-")
        mapped_start = mapping[int(start)]
        mapped_end = mapping[int(end)]
        processed_ranges.append([mapped_start, mapped_end])
    
    return processed_ranges

# def mask_ranges(aes_range, masks):
#     range_low, range_high = aes_range
#     aes_range_list = range(range_low, range_high)
    
#     new_aes_range = []
#     for [mask_low, mask_high] in masks:
#         if mask_low in aes_range_list and mask_high in aes_range_list:
#             new_aes_range.append([range_low, mask_low - 1])
#             new_aes_range.append([mask_high + 1, range_high])
            
#         else:
#             if mask_low in aes_range_list and mask_low < range_high:
#                 range_high = mask_low - 1
                
#             if mask_high in aes_range_list and mask_high > range_low:
#                 range_low = mask_high + 1
            
#             new_aes_range.append([range_low, range_high])
            
#     return new_aes_range
    
    
def mask_ranges(aes_range, masks):
    range_low, range_high = aes_range
    new_aes_range = [[range_low, range_high]]  # Start with the entire range

    for mask_low, mask_high in masks:
        updated_aes_range = []
        
        for r_low, r_high in new_aes_range:
            # No overlap with the current range
            if mask_high < r_low or mask_low > r_high:
                updated_aes_range.append([r_low, r_high])
            else:
                # Check the part before the mask
                if mask_low > r_low:
                    updated_aes_range.append([r_low, mask_low - 1])
                # Check the part after the mask
                if mask_high < r_high:
                    updated_aes_range.append([mask_high + 1, r_high])
        
        new_aes_range = updated_aes_range  # Update the list of ranges
    
    return new_aes_range

def rotate_aes(aes_list, order):
    
    for _ in range(order):
        last = aes_list.pop()
        aes_list.insert(0, last)
        
    return aes_list

def process_aes_df(aes_csv_path, mapping, augment=False):
    aes_df = pd.read_csv(aes_csv_path, header=None, names=["aes", "middle", "ranges"])
    
    if "mask" in aes_df["aes"].tolist():
        msa_map = {i:i for i in range(10**4)} # dirty patch
        masked_ranges = aes_df[aes_df["aes"] == "mask"]["ranges"].apply(lambda x: process_range(x, msa_map)).tolist()[0]
    else:
        masked_ranges = []
    
    aes_df = aes_df[aes_df["aes"] != "mask"]
    aes_df["ranges"] = aes_df["ranges"].apply(lambda x: process_range(x, mapping))
    
    if len(masked_ranges) > 0:
        aes_df_exploded = aes_df[aes_df["aes"] != "mask"].explode("ranges")
        
        aes_df_exploded["ranges"] = aes_df_exploded["ranges"].apply(lambda x: mask_ranges(x, masked_ranges))
        aes_df = aes_df_exploded.groupby("aes").agg({"ranges": 'sum'}).reset_index()
    
    aes_mapping = aes_df.set_index("aes")["ranges"].to_dict()
    
    if augment:
        augmented_mapping = {}
        for k, v in aes_mapping.items(): 
            augmented_mapping[f"{k}_0"] = v
            for order in range(1, len(v)):
                augmented_mapping[f"{k}_{order}"] = rotate_aes(v[:], order)
        
        aes_mapping = augmented_mapping
        
    return aes_mapping


def remove_inter_aes_bonds(dot_bracket):
    stack = []
    illegal_indexes = []
    for index, symbol in enumerate(dot_bracket):
        if symbol == "(":
            stack.append(index)
        elif symbol == ")":
            if len(stack) > 0:
                stack.pop()
            else:
                illegal_indexes.append(index)
    
    illegal_indexes.extend(stack)
    dot_bracket = list(dot_bracket)
    
    for i in illegal_indexes:
        dot_bracket[i] = '.'
    
    return "".join(dot_bracket)

def create_dot_bracket(bps, max_len):
    dot_bracket = ["."] * max_len
    for start, stop in bps:
        start, stop = sorted([start, stop])
        dot_bracket[start] = "("
        dot_bracket[stop] = ")"
    
    return "".join(dot_bracket)

def break_msa(fasta_path, aes_mapping, aes_bp_mapping):
    aes_records = {}
    
    for aes_name, ranges in aes_mapping.items():
        modified_records = []
        for record in SeqIO.parse(fasta_path, "fasta"):
            new_id = record.description.replace(" ", "__")
            sequence_records = ""
            current_dot_bracket = ""
            for [start, stop] in ranges:
                sequence_records += record.seq[int(start)-1:int(stop)]
                # current_dot_bracket += dot_bracket[int(start):int(stop)+1]
                
            new_record = SeqRecord(sequence_records, id=new_id, description="")
            modified_records.append(new_record)
        
        current_dot_bracket = create_dot_bracket(aes_bp_mapping[aes_name], len(sequence_records))
        # current_dot_bracket = remove_inter_aes_bonds(current_dot_bracket)
        
        aes_records[aes_name] = [modified_records, current_dot_bracket]
    
    return aes_records

def write_stockholm(aes_records, base_path):
    stockholm_paths = []
    for aes, [records, dot_bracket] in aes_records.items():
        sto_path = os.path.join(base_path, f"aes_{aes}.stockholm")
        os.makedirs(os.path.dirname(sto_path), exist_ok=True)
        SeqIO.write(records, sto_path, "stockholm")
        
        insert_line_to_file(sto_path, f"#=GC SS_cons {dot_bracket}", 1)
        stockholm_paths.append(sto_path)
        
    return stockholm_paths

def write_fasta(aes_records, base_path):
    fasta_paths = []
    for aes, [records, dot_bracket] in aes_records.items():
        sto_path = os.path.join(base_path, f"aes_{aes}", "msa.fasta")
        os.makedirs(os.path.dirname(sto_path), exist_ok=True)
        SeqIO.write(records, sto_path, "fasta")
        fasta_paths.append(sto_path)
        
    return fasta_paths

def insert_line_to_file(file_path, line, index):
    with open(file_path, "r") as f:
        contents = f.readlines()
        
    contents.insert(index, line + "\n")
    
    with open(file_path, "w") as f:
        contents = "".join(contents)
        f.write(contents)
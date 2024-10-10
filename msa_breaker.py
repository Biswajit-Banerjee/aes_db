from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from utils import log, modify_extension
import pandas as pd
from pathlib import Path
import os
import json

def process_range(range_str, mapping):
    processed_ranges = []
    range_pairs = range_str.split(";")
    
    for pair in range_pairs:
        start, end = pair.split("-")
        mapped_start = mapping[int(start)]
        mapped_end = mapping[int(end)]
        processed_ranges.append([mapped_start, mapped_end])
    
    return processed_ranges

def process_aes_df(aes_csv_path, mapping):
    # Read the CSV file
    aes_df = pd.read_csv(aes_csv_path, header=None, names=["aes", "middle", "ranges"])
    
    aes_df["ranges"] = aes_df["ranges"].apply(lambda x: process_range(x, mapping))
    
    # Return the result as a dictionary
    return aes_df.set_index("aes")["ranges"].to_dict()

# def process_aes_df(aes_csv_path, mapping_path):
#     aes_df = pd.read_csv(aes_csv_path, header=None)
#     aes_df.columns = ["aes", "middle", 'ranges']
    
#     with open(mapping_path, 'r') as fp:
#         mapping = json.load(fp)
    
#     aes_df.ranges = aes_df.ranges.apply(lambda x: list(map(lambda y: [mapping[y.split("-")[0]], mapping[y.split("-")[1]]], x.split(";"))))    
#     aes_df.set_index("aes", inplace=True)
    
#     return aes_df["ranges"].to_dict() 


def remove_inter_aes_bonds(dot_bracket):
    stack = []
    
    illigal_indexes = []
    for index, symbol in enumerate(dot_bracket):
        if symbol == "(":
            stack.append(index)
            
        elif symbol == ")":
            if len(stack) > 0:
                stack.pop()
            else:
                illigal_indexes.append(index)
    
    illigal_indexes.extend(stack)
    dot_bracket = list(dot_bracket)
    
    for i in illigal_indexes:
        dot_bracket[i] = '.'
    
    return "".join(dot_bracket)
        

def break_msa(fasta_path, aes_mapping, dot_bracket = ""):

    aes_records = {}
    # lines = []
    
    for aes_name, ranges in aes_mapping.items():
        
        modified_records = []
        
        for record in SeqIO.parse(fasta_path, "fasta"):
            new_id = record.description.replace(" ", "__")
            sequence_records = ""
            current_dot_bracket = ""
            for [start, stop] in ranges:
                sequence_records += record.seq[int(start)-1:int(stop)]
                current_dot_bracket += dot_bracket[int(start):int(stop)+1]
                
            new_record = SeqRecord(sequence_records, id=new_id, description="")
            modified_records.append(new_record)
            
            # if "Escherichia" in new_id:
            # lines.extend([str(record.description) + f" AES{aes_name}", current_dot_bracket, str(sequence_records), "\n"])
        
            # lines.extend([f"'{str(record.description)} AES{aes_name}': '{str(sequence_records)}',"])    
        
        # fix the dot bracket
        current_dot_bracket = remove_inter_aes_bonds(current_dot_bracket)
        aes_records[aes_name] = [modified_records, current_dot_bracket]
        # lines.append(current_dot_bracket)
    
    # with open("test_file.txt", "w") as file:
    #     file.writelines("\n".join(lines))
        
    
    return aes_records

def write_stockholm(aes_records):
    base_path = Path(os.environ.get("CM_BASE", "./results"))
    sotckholm_paths = []
    for aes, [records, dot_bracket] in aes_records.items():
        sto_path = os.path.join(base_path, f"aes_{aes}", "msa.stockholm")
        os.makedirs(os.path.dirname(sto_path), exist_ok=True)
        SeqIO.write(records, sto_path, "stockholm")
        
        # add secondary structure to the file
        insert_line_to_file(sto_path, f"#=GC SS_cons {dot_bracket}", 1)
        sotckholm_paths.append(sto_path)
        
    return sotckholm_paths
    

def insert_line_to_file(file_path, line, index):
    with open(file_path, "r") as f:
        contents = f.readlines()
        
    contents.insert(index, line + "\n")
    
    with open(file_path, "w") as f:
        contents = "".join(contents)
        f.write(contents)

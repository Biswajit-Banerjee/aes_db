import subprocess
import os
from .utils import log, modify_extension
from .config import Binaries as bins

def run_cmbuild(sto_path, postfix=""):
    """
    Run the cmbuild tool on a stockholm file.
    
    cmbuild -F <cm output file path> <stockholm input file path>

    Args:
        sto_path (str): The path to the Stockholm file.
        inf_path (str, optional): Path to infernal. 
        postfix (str, optional): Postfix to append to the output directory name. Defaults to "".

    Returns:
        str: The path to the generated cmbuild file if successful, None otherwise.
    """
    cm_path = modify_extension(sto_path, "cm").replace("stockholms", "covariance_models" + postfix)    
    os.makedirs(os.path.split(cm_path)[0], exist_ok=True)
    result = subprocess.run([bins.INF_115 / "src/cmbuild", "-F", cm_path, sto_path], 
                            capture_output=True, text=True)
    
    if result.returncode != 0:
        log(f"Error running cmbuild for {sto_path}, see below for more details", "error")
        print(result.stderr)
        return None
    
    return cm_path

def run_cm_convert(cm_path, flag="--mlhmm"):
    """
    Run the hmmbuild tool on a stockholm file.
    
    hmmbuild -F <hmm output file path> <stockholm input file path>

    Args:
        sto_path (str): The path to the Stockholm file.
    Returns:
        str: The path to the generated hmmbuild file if successful, None otherwise.
    """
    if flag == "--mlhmm":
        file_path = modify_extension(cm_path, "hmm").replace("covariance_models", "hmm_profiles")  
    elif flag == "-1":
        file_path = modify_extension(cm_path, "cm").replace("covariance_models", "cm_1p0p2")  
    else:
        raise Exception("Not implemented Yet!")
    
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    
    command = [str(bins.CM_CONVERT), flag, cm_path]
    result = subprocess.run(
        command, 
        capture_output=True, text=True
    )
    
    if result.returncode != 0:
        log(f"Error running hmmbuild for {cm_path}, see below for more details", "error")
        print(result, "\n", " ".join(command))
        return None
    else:
        with open(file_path, 'w') as f:
            f.writelines(result.stdout)
    return file_path
    

def extract_base_pairs(bp_json, index_mapping=None):
    annotations = bp_json["annotations"]
    
    base_pairs = []
    marked_neucleotides = []
    
    for value in annotations:                
        if value["bp"] == "cWW" and value["crossing"] == "0":
            start, stop = int(value['seq_id1']), int(value['seq_id2'])
            if index_mapping:
                start, stop = index_mapping[start], index_mapping[stop]
            
            if start not in marked_neucleotides and stop not in marked_neucleotides:   
                base_pairs.append([start, stop])
                marked_neucleotides.extend([start, stop])
            else:
                print([start, stop])
    
    return base_pairs

def continous_aes_mapping(aes_mapping):
    continous_map = {}
    for aes, ranges in aes_mapping.items():
        temp_list = []
        for r in ranges:
            temp_list.extend(list(range(r[0], r[1]+1)))
        
        continous_map[aes] = temp_list
    
    return continous_map

def colocate_basepairs_in_aes(base_pairs, aes_mapping):
    aes_bp_mapping = {}
    continous_map = continous_aes_mapping(aes_mapping)
    
    temp_bp = base_pairs[:]
    for aes, indexes in continous_map.items():
        aes_bp_list = []
        print(f"AES: {aes} Available BPs: {len(temp_bp)}")
        for start, stop in temp_bp[:]:
            if start in indexes and stop in indexes:
                local_start = indexes.index(start)
                local_stop  = indexes.index(stop)
                aes_bp_list.append([local_start, local_stop])
                # temp_bp.remove([start, stop])
        
        aes_bp_mapping[aes] = aes_bp_list
    print(f"Couldn't assign BPs: {len(temp_bp)}")
    
    return aes_bp_mapping    

    
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
    mapping = {}
    
    current_index = start_index
    for i, symbol in enumerate(gapped_sequence, start_index):
        if symbol != "-":
            mapping[current_index] = i
            current_index += 1
    
    return mapping
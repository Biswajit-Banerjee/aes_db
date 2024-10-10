import subprocess 
from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
import multiprocessing as mp
from .utils import modify_extension

def fasta_to_clustral(fasta_path, output_file):
    
    if not output_file:
        output_file = modify_extension(fasta_path, "aln")
    
    # Read the FASTA alignment
    alignment = AlignIO.read(fasta_path, "fasta")

    # Write the alignment in CLUSTAL format
    with open(output_file, "w") as f:
        AlignIO.write(alignment, f, "clustal")
    
    return output_file

def run_ipknot(fasta_path, output_file=""):
    output_file = fasta_to_clustral(fasta_path, output_file)
    n_threads = mp.cpu_count() - 1
    
    output = subprocess.run(["/home/sumon/repos/aes_db/ipknot", "-n", str(n_threads), output_file], capture_output=True, text=True)
    _, probable_seq, dot_bracket = output.stdout.strip().split("\n")
    return probable_seq, dot_bracket

if __name__ == "__main__":
    fasta = "/home/sumon/repos/aes_db/cm_builder_storage/results/aes_30/msa.fasta"
    probable_seq, dot_bracket = run_ipknot(fasta)

    


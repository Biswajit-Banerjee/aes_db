from rna_cm_builder.core import process_single_msa, process_single_msa_using_ipknot


# msa_data = {
#     "fasta_path": "/home/sumon/repos/aes_db/data/all_fas/bactRNAseP.fas",
#     "organism": "Structure sequence", # structure seq
#     "pdb_id": "3Q1Q", # 3Q1Q
#     "chain": "B", # B
#     "aes_csv_path": "/home/sumon/repos/aes_db/data/aes_defs/TM_RnasP_AES_Defs.csv" # create 
#   }

# process_single_msa(msa_data)

msa_data = {
    "fasta_path": "/home/sumon/repos/aes_db/data/msa/23S.fas",
    "organism": "Escherichia coli",
    "pdb_id": "4V9D", # 3Q1Q
    "chain": "CA", # B
    "aes_csv_path": "/home/sumon/repos/aes_db/data/aes_defs/AES_defs_LSUb.csv", # create 
    "out_dir": "./test"
  }

process_single_msa(msa_data)

# msa_data = {
#     "fasta_path": "/home/sumon/repos/aes_db/data/all_fas/16S.fas",
#     "organism": "Escherichia coli",
#     "pdb_id": "4V9D", # 3Q1Q
#     "chain": "AA", # B
#     "aes_csv_path": "/home/sumon/repos/aes_db/data/aes_defs/aes_defs_ssu.csv", # create 
#     "out_dir": "./test"
#   }

# process_single_msa(msa_data)

from Bio import AlignIO
from Bio.AlignIO.PhylipIO import SequentialPhylipWriter

data_path = "../data/Wheat_Relative_History_Data_Glemin_et_al/"
concat_file = "triticeae_allindividuals_OneCopyGenes.fasta"

fasta_file = data_path + concat_file
phylip_file = "../results/10-triticeae_allindividuals_OneCopyGenes.phylip"

# Load the alignment
with open(fasta_file, "r") as f_in:
    alignment = AlignIO.read(f_in, "fasta")

# Find the length of the longest sequence ID to prevent truncation
max_id_len = max(len(record.id) for record in alignment)

# Write out the alignment
with open(phylip_file, "w") as f_out:
    writer = SequentialPhylipWriter(f_out)
    writer.write_alignment(alignment, id_width=max_id_len + 3) # 3 additional padding


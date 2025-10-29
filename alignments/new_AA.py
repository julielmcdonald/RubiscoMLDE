pip install biopython pandas

import re
from Bio import AlignIO
import pandas as pd

# === INPUT ===
#replace with your aligned fasta file
alignment_file = "/Users/juliemcdonald/Documents/Shoulders_lab/Rubisco_data/AA_identities/rhodophyta_rbcL_aligned.fasta"

#replace with AA site you want to test (should correspond to the site in the first sequence of your alignment)
unaligned_position = 228

alignment = AlignIO.read(alignment_file, "fasta")

# === STEP 1: Get aligned index of position in first sequence ===
first_seq = alignment[0].seq
non_gap_counter = 0
aligned_index = None

for i, aa in enumerate(first_seq):
    if aa != "-":
        non_gap_counter += 1
    if non_gap_counter == unaligned_position:
        aligned_index = i
        break

if aligned_index is None:
    raise ValueError("Could not find position in first sequence.")

print(f"Aligned column corresponding to position: {aligned_index}")

# === STEP 2: Extract amino acid & species name from all sequences ===
data = []

for record in alignment:
    seq_id = record.id
    description = record.description
    aa = record.seq[aligned_index]

    # Try to extract species name in brackets
    match = re.search(r"\[([^\]]+)\]", description)
    species = match.group(1) if match else "Unknown"

    data.append({
        "Sequence_ID": seq_id,
        "Species": species,
        "AminoAcid_at_pos": str(aa)
    })
    

# === STEP 3: Convert to DataFrame ===
df = pd.DataFrame(data)

aa_counts = df["AminoAcid_at_pos"].value_counts().reset_index()
aa_counts.columns = ["AminoAcid", "Count"]
print(aa_counts)

#Determine which species have a given amino acid at given position

#change to desired AA
input_AA = "R"

AA_seqs = df[df["AminoAcid_at_pos"] == input_AA]
print("\nSequences with",input_AA,"at position",unaligned_position)
print(AA_seqs[["Sequence_ID", "Species"]])
table = AA_seqs[["Sequence_ID", "Species"]]

# === Save results ===
#df.to_csv(" ", index=False)
#aa_counts.to_csv(" ", index=False)
#table.to_csv(" ",index=False)

import pandas as pd

df = pd.read_csv("list_cps_genes_proteins_sequences.csv")

#loop through each unique gene
for gene in df['A'].unique():
    #subset sequences for this gene
    sequences = df[df['A'] == gene]['B']
    
    #create .faa file for the gene
    with open(f"{gene}.faa", "w") as f:
        for i, seq in enumerate(sequences, start=1):
            # Add stop codon '*' if not already present
            if not seq.endswith("*"):
                seq += "*"
            f.write(f">{gene}_{i}\n{seq}\n")

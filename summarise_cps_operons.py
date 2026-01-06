import os
from collections import defaultdict
from Bio import SeqIO
import pandas as pd

def summarize_genbank_features(folder):
    gene_data = defaultdict(lambda: {
        'files': set(),
        'count': 0,
        'locus_tags': set(),
        'notes': set(),
        'longest_translation': ''
    })

    files = [f for f in os.listdir(folder) if f.lower().endswith(('.gb', '.gbk'))]

    for file in files:
        filepath = os.path.join(folder, file)
        seen_in_file = set()

        for record in SeqIO.parse(filepath, "genbank"):
            for feature in record.features:
                if feature.type in ['gene', 'CDS']:
                    if 'gene' in feature.qualifiers:
                        gene_name = feature.qualifiers['gene'][0]
                    elif 'locus_tag' in feature.qualifiers:
                        gene_name = feature.qualifiers['locus_tag'][0]
                    elif 'product' in feature.qualifiers:
                        gene_name = feature.qualifiers['product'][0]
                    else:
                        continue
                    
                    gene_name = gene_name.strip().upper()

                    if (gene_name, file) not in seen_in_file:
                        gene_data[gene_name]['count'] += 1
                        seen_in_file.add((gene_name, file))

                    gene_data[gene_name]['files'].add(file)

                    if 'locus_tag' in feature.qualifiers:
                        gene_data[gene_name]['locus_tags'].update(feature.qualifiers['locus_tag'])

                    if 'note' in feature.qualifiers:
                        gene_data[gene_name]['notes'].update(feature.qualifiers['note'])

                    if 'translation' in feature.qualifiers:
                        translation = feature.qualifiers['translation'][0]
                        if len(translation) > len(gene_data[gene_name]['longest_translation']):
                            gene_data[gene_name]['longest_translation'] = translation

    gene_summary = []
    for gene, info in gene_data.items():
        gene_summary.append({
            'Gene': gene,
            'Files': ';'.join(sorted(info['files'])),
            'Files_count': len(info['files']),
            'Total_count': info['count'],
            'Locus_tags': ';'.join(sorted(map(str, info['locus_tags']))),
            'Longest_translation': info['longest_translation'],
            'Notes': ' || '.join(map(str, info['notes'])),
        })

    df = pd.DataFrame(gene_summary)
    df.to_csv('gene_summary_single_longest_translation.csv', index=False)
    print("Summary saved as 'gene_summary_single_longest_translation.csv'")

folder_path = '/x/pneumococcal_cps-reference_genebank'
summarize_genbank_features(folder_path)

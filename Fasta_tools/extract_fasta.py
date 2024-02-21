import argparse
import re
from Bio import SeqIO

def filter_and_write_fasta(input_file, output_file, gene_names_file, longest):
    # Read the list of gene names
    with open(gene_names_file, 'r') as file:
        gene_names = [line.strip() for line in file]

    # Store found sequences ({gene_name: sequence_object})
    found_sequences = {}

    # Read and filter the FASTA file
    for record in SeqIO.parse(input_file, "fasta"):
        for gene_name in gene_names:
            # Use regular expression for fuzzy matching, allowing extra characters (e.g., .1, .2) after the gene name
            if re.match(f"{gene_name}(\.\d+)?$", record.id):
                if gene_name not in found_sequences:
                    found_sequences[gene_name] = [record]
                else:
                    found_sequences[gene_name].append(record)

    # If --longest is specified, keep only the longest sequence for each gene match
    if longest:
        for gene_name in found_sequences:
            found_sequences[gene_name] = [max(found_sequences[gene_name], key=lambda x: len(x.seq))]

    # Write the output file
    with open(output_file, "w") as output_handle:
        for gene_list in found_sequences.values():
            for record in gene_list:
                SeqIO.write(record, output_handle, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter FASTA file for specific genes and optionally get the longest sequence.")
    parser.add_argument("input_file", help="Input FASTA file path")
    parser.add_argument("gene_names_file", help="File containing gene names to filter by")
    parser.add_argument("output_file", help="Output file path for filtered FASTA")
    parser.add_argument("--longest", action="store_true", help="Only output the longest sequence for each gene match")

    args = parser.parse_args()

    filter_and_write_fasta(args.input_file, args.output_file, args.gene_names_file, args.longest)


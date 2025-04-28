import json
from Bio import SeqIO
import pandas as pd
import click

def create_chain_map(fasta_file, json_file, output_file):
    """Create chain mapping TSV from FASTA and AlphaFold3 JSON"""
    try:
        # Parse FASTA headers
        fasta_entries = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            header = record.description
            if "#" in header:
                name_part, copies_part = header.split("#", 1)
                copies = int(copies_part.split()[0])
                name = name_part.split()[0].strip()
            else:
                name = header.split()[0].strip()
                copies = 1
            fasta_entries.append((name, copies))

        # Parse JSON chains
        with open(json_file) as f:
            af_data = json.load(f)

        json_chains = []
        for seq in af_data["sequences"]:
            for seq_type in ["protein", "dna", "rna", "ligand"]:
                if seq_type in seq:
                    json_chains.extend(seq[seq_type]["id"])
                    break

        # Validate chain counts
        total_fasta_chains = sum(copies for _, copies in fasta_entries)
        if total_fasta_chains != len(json_chains):
            raise ValueError(
                f"Chain count mismatch: FASTA expects {total_fasta_chains} chains, "
                f"JSON contains {len(json_chains)} chains"
            )

        # Create mapping
        mapping = []
        chain_idx = 0
        for name, copies in fasta_entries:
            for copy_num in range(1, copies + 1):
                mapping.append({
                    "ChainID": json_chains[chain_idx],
                    "ChainName": name,
                    "CopyNumber": copy_num,
                    "TotalCopies": copies
                })
                chain_idx += 1

        # Save TSV
        pd.DataFrame(mapping).to_csv(output_file, sep="\t", index=False)
        return True

    except Exception as e:
        raise RuntimeError(f"Mapping error: {str(e)}")

@click.command()
@click.argument("fasta_file", type=click.Path(exists=True))
@click.argument("json_file", type=click.Path(exists=True))
@click.option("-o", "--output", default="chain_mapping.tsv", help="Output TSV file")
def map_chains(fasta_file, json_file, output):
    """
    Generate chain-to-sample mapping for AlphaFold3 results.

    Arguments:
    - fasta_file: Input FASTA file created by af3build.
                  The FASTA file must include headers with the format:
                  >Name#Copies [Modifications] [OriginalHeader]
    - json_file: Input AlphaFold3 JSON file created by af3tools fasta2json.
                 The JSON file must include sequences with assigned chain IDs.
    - output: Name of the output TSV file (default: chain_mapping.tsv).
    """
    try:
        if create_chain_map(fasta_file, json_file, output):
            click.secho(f"Successfully created chain mapping: {output}", fg="green")
    except Exception as e:
        click.secho(str(e), fg="red")
        raise click.Abort()

if __name__ == "__main__":
    map_chains()

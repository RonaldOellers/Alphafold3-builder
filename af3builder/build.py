import click
from pathlib import Path
from .core import AF3Builder
from .exceptions import AF3Error

@click.command()
@click.argument("input_file", type=click.Path(exists=True), help="Provide input yaml or tsv file")
@click.option("-em", "--email", help="Email for NCBI access")
@click.option("-o", "--output", help="Output FASTA file", default="af3_input.fasta")
@click.option("-v", "--verbose", is_flag=True, help="Verbose output for debugging")

def build(input_file, email, output, verbose):
    """Convert biological data to AlphaFold3 input FASTA"""
    try:
        builder = AF3Builder(email=email)
        builder.build(input_file, output, verbose)
        click.secho(f"Successfully created {output}", fg="green")
    except AF3Error as e:
        click.secho(f"Error: {str(e)}", fg="red")
        raise click.Abort()

if __name__ == "__main__":
    build()

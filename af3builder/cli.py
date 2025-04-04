import click
from pathlib import Path
from .core import AF3Builder
from .exceptions import AF3Error

@click.command()
@click.argument("input_file", type=click.Path(exists=True))
@click.option("--email", help="Email for NCBI access")
@click.option("-o", "--output", default="af3_input.fasta")
def main(input_file, email, output):
    """Convert biological data to AlphaFold3 input"""
    try:
        builder = AF3Builder(email=email)
        builder.build(input_file, output)
        click.echo(f"Successfully created {output}")
    except AF3Error as e:
        click.secho(f"Error: {str(e)}", fg="red")
        raise click.Abort()

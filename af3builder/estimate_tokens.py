import click
from .core import AlphaFold3TokenCounter
from .exceptions import AF3Error

@click.command()
@click.argument('fasta_path', type=click.Path(exists=True))
@click.option('--verbose', is_flag=True, help='Verbose output')
@click.option('--recommendedgpu', is_flag=True, help='Show GPU recommendations')
@click.option('--smile_leniency', is_flag=True, 
              help='Add 5% token headroom for ligand/small molecule complexity')
def estimateTokens(fasta_path, verbose, recommendedgpu, smile_leniency):
    """AlphaFold3 token estimator with hardware recommendations"""
    try:
        counter = AlphaFold3TokenCounter(fasta_path, verbose, recommendedgpu, smile_leniency)
        counter.parse()
        counter.summary()
        click.secho("Token estimation completed successfully.", fg="green")
    except AF3Error as e:
        click.secho(f"Error: {str(e)}", fg="red")
        raise click.Abort()
    except Exception as e:
        click.secho(f"Unexpected error: {str(e)}", fg="red")
        raise click.Abort()

if __name__ == "__main__":
    estimateTokens()

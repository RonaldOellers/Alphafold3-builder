# AlphaFold3 Input Preparation Toolkit

This toolkit converts biological data (TSV/YAML) into AlphaFold3-compatible FASTA files. It supports raw sequences, database IDs (UniProt/NCBI), ligands, SMILES strings, and post-translational modifications.

## Features
- **Input Formats**: TSV/YAML with fields for ID, Type, Copies, Modifications.
- **Database Support**: UniProt (proteins), NCBI (DNA/RNA).
- **FASTA Output**:
  - Proper formatting for AlphaFold3 tools.
  - Automatic DNA→RNA conversion.
  - Post-translational modifications (`&position_code`).
  - Oligomer counts (`#copies`).
  - Separation between database ID info and custom FASTA headers.

## Setup
### Using Conda Environment
```bash
conda env create -f envs/af3env.yaml
conda activate af3env
pip install -e .
```

### Using Pytest Environment for Testing
```bash
conda env create -f envs/pytest_env.yaml
conda activate af3_pytest
pytest tests/
```


## Usage
### Input File Format (TSV Example)
```TSV
ID Type Copies Modifications Name
P0DTD1 protein 2 &199_CSO SpikeProtein
ligand#3 ligand 1 MG LigandMG
NM_001301244 rna 1 &199_PSU MyRNA
ATGGTACCTA dna 1 CustomDNA CustomGeneName
```

### Input File Format (YAML Example)
```yaml
- ID: P0DTD1
  TYPE: protein
  COPIES: 2
  MODIFICATIONS: &199_CSO
  NAME: SpikeProtein

- ID: ligand
  TYPE: ligand
  COPIES: 1
  MODIFICATIONS: MG
  NAME: LigandMG

- ID: NM_001301244
  TYPE: rna
  COPIES: 1
  MODIFICATIONS: &199_PSU

- ID: ATGGTACCTA
  TYPE: dna
  COPIES: 1
  NAME: CustomGeneName

```

### Extended Format explanation
```yaml
# Database-derived protein with modification
- ID: P0DTD1               # UniProt ID
  TYPE: protein
  COPIES: 2                # Homodimer
  MODIFICATIONS: &199_CSO  # Cysteine sulfenic acid at position 199
  NAME: SpikeProtein       # Custom display name

# Ligand entry (raw input)
- ID: ligand
  TYPE: ligand
  COPIES: 1
  MODIFICATIONS: MG        # Magnesium ion
  NAME: LigandMG           # Required for ligands

# Database-derived RNA with modification
- ID: NM_001301244         # NCBI accession
  TYPE: rna
  MODIFICATIONS: &199_PSU  # Pseudouridine at position 199
  # Name omitted → uses ID in FASTA header

# Custom DNA sequence
- ID: ATGGTACCTA           # Raw sequence
  TYPE: dna
  NAME: CustomGeneName     # Required for custom sequences
```

### Command to Generate FASTA File
```bash
af3build --email your@lab.com -o af3_input.fasta input.tsv
```

### Output FASTA Example (with separation between DB info and custom headers)
```FASTA
>P0DTD1
>proteinSpikeProtein &199_CSO #2
MAEGASTERDA...

>ligandLigandMG
MG

>NM_001301244
>rnaMyRNA &199_PSU
AUGGCUAAU...

>dnaCustomGeneName
ATGGTACCTA...
```

## Input Specifications

| Column         | Required | Values                |
|----------------|----------|-----------------------|
| `ID`           | Yes      | Database ID or raw sequence |
| `TYPE`         | Yes      | `protein`, `dna`, `rna`, `ligand`, `smile` |
| `Copies`       | No       | Integer (default=1)   |
| `Modifications`| No       | `&position_code`      |
| `Name`         | Required for custom sequences or ligands | String |

## Entrez Email Requirement

For fetching sequences from NCBI databases (DNA/RNA), you need an Entrez email. You can sign up for an account at [NCBI](https://www.ncbi.nlm.nih.gov/account/) and use your registered email address when running the program.

## Related Tools

- [AlphaFold3_tools](https://github.com/snufoodbiochem/Alphafold3_tools): Includes `fasta2json.py` for converting FASTA to JSON.
- [AlphaFold3 Official](https://github.com/google-deepmind/alphafold3): Main repository for AlphaFold3.


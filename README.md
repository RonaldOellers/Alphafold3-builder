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
ID	TYPE	COPIES	MODIFICATIONS	NAME
P0DTD1	protein	2	&199_CSO	SpikeProtein
ligand	ligand	1	MG	LigandMG
NM_001301244	rna	1	&199_PSU	NM_001301244
ATGGTACCTA	dna	1		CustomGeneName
```

### Input File Format (YAML Example)
```yaml
- ID: P0DTD1
  TYPE: protein
  COPIES: 2
  MODIFICATIONS: &199_CSO
  NAME: Covid_SpikeProtein

- ID: MG
  TYPE: ligand
  COPIES: 1
  MODIFICATIONS:
  NAME: Magnesium

- ID: NM_001301244
  TYPE: rna
  COPIES: 3
  MODIFICATIONS: &199_PSU
  Name: NM_001301244

- ID: ATGGTACCTA
  TYPE: dna
  COPIES: 2
  NAME: CustomGeneName

```

### Extended Format explanation
```yaml
# Database-derived protein with modification
- ID: P0DTD1               # UniProt ID
  TYPE: protein
  COPIES: 2                # Homodimer
  MODIFICATIONS: &199_CSO  # Cysteine sulfenic acid at position 199
  NAME: SpikeProtein       # Custom display name (Always recomended)

# Ligand entry (raw input)
- ID: ligand
  TYPE: ligand
  COPIES: 1
  MODIFICATIONS: MG        # Magnesium ion
  NAME: LigandMG

# Database-derived RNA with modification
- ID: NM_001301244         # NCBI accession
  TYPE: rna
  COPIES: 3
  MODIFICATIONS: &199_PSU  # Pseudouridine at position 199
  # Name omitted → leaves empty position in FASTA header

# Custom DNA sequence
- ID: ATGGTACCTA           # Raw sequence
  TYPE: dna
  COPIES: 2
  NAME: CustomGeneName
```

### Command to Generate FASTA File
```bash
af3build --email your@lab.com -o af3_input.fasta input.tsv
```

### Output FASTA Example (with separation between DB info and custom headers)
```FASTA
>Covid_SpikeProtein :: #2 :: protein :: &199_CSO :: sp|P0DTD1|R1AB_SARS2 Replicase polyprotein 1ab OS=Severe acute respiratory syndrome coronavirus 2 OX=2697049 GN=rep PE=1 SV=1
MESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVE...

>Magnesium :: ligand
MG

>NM_001301244 :: #3 :: rna :: &199_PSU :: NM_001301244.2 Homo sapiens tropomyosin 1 (TPM1), transcript variant Tpm1.2, mRNA
GCUCGCACUCCCGCUCCUCCGCCCGACCGCGCGCUCG...

>CustomGeneName :: #2 :: dna
ATGGTACCTA
```

## Input Specifications

| Column         | Required | Values                |
|----------------|----------|-----------------------|
| `ID`           | Yes      | Database ID or raw sequence |
| `TYPE`         | Yes      | `protein`, `dna`, `rna`, `ligand`, `smile` |
| `Copies`       | No       | Integer (default=1)   |
| `Modifications`| No       | `&position_code`      |
| `Name`         | Recommended | String |

## Entrez Email Requirement

For fetching sequences from NCBI databases (DNA/RNA), you need an Entrez email. You can sign up for an account at [NCBI](https://www.ncbi.nlm.nih.gov/account/) and use your registered email address when running the program.

## Related Tools

- [AlphaFold3_tools](https://github.com/snufoodbiochem/Alphafold3_tools): Includes `fasta2json.py` for converting FASTA to JSON.
- [AlphaFold3 Official](https://github.com/google-deepmind/alphafold3): Main repository for AlphaFold3.


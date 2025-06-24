# AlphaFold3 Input Preparation Toolkit
[![AF3Builder CI](https://github.com/RonaldOellers/Alphafold3-builder/actions/workflows/CI.yaml/badge.svg?branch=main&event=push)](https://github.com/RonaldOellers/Alphafold3-builder/actions/workflows/CI.yaml)

Create Alphafold3 input using FASTA and prepare conversion to JSON.

This toolkit converts sample sheets (TSV/YAML) into [AlphaFold3_tools](https://github.com/snufoodbiochem/Alphafold3_tools)-compatible FASTA files. It supports raw sequences, database IDs (UniProt/NCBI), ligands, SMILES strings, and post-translational modifications.

These FASTA files can easily be converted to the JSON format required by AlphaFold3 using [AlphaFold3_tools](https://github.com/snufoodbiochem/Alphafold3_tools).

It also maps the Entries in the FASTA and JSON to the Chains as they will appear in the AlphaFold3 Output.

## Features
- **Input Formats**: TSV/YAML with fields for ID, Type, Copies, Modifications.
- **Database Support**: UniProt (proteins), NCBI (DNA/RNA).
- **FASTA Output**:
  - Proper formatting for AlphaFold3 tools.
  - Automatic DNA→RNA conversion.
  - Post-translational modifications (`&position_code`).
  - Oligomer counts (`#copies`).
  - Separation between database ID info and custom FASTA headers.
- **Chains Map**: Get the Order of Chains with names and settings for reference in downstream analysis

## Setup
### Using Conda Environment
```bash
conda env create -f envs/af3env.yaml #or any env with python 3.12
conda activate AF3builder
pip install -e .
```

### Using Pytest Environment for Testing (unfinished)
```bash
conda env create -f envs/pytest_env.yaml
conda activate af3_pytest
pytest tests/
```


# Usage
## af3build
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
  NAME: SpikeProtein       # Custom display name (Always recommended)

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
  # NAME omitted → leaves an empty position in FASTA header

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
>Covid_SpikeProtein &199_CSO sp|P0DTD1|R1AB_SARS2 Replicase polyprotein 1ab OS=Severe acute respiratory syndrome coronavirus 2 OX=2697049 GN=rep PE=1 SV=1 #2
MESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVE...

>Magnesium
MG

>NM_001301244 &199_PSU NM_001301244.2 Homo sapiens tropomyosin 1 (TPM1), transcript variant Tpm1.2, mRNA #3
GCUCGCACUCCCGCUCCUCCGCCCGACCGCGCGC...

>CustomGeneName #2
ATGGTACCTA
```

### Input Specifications

| Column         | Required | Values                |
|----------------|----------|-----------------------|
| `ID`           | Yes      | Database ID or raw sequence |
| `TYPE`         | Yes      | `protein`, `dna`, `rna`, `ligand`, `smile` |
| `Copies`       | No       | Integer (default=1)   |
| `Modifications`| No       | `&position_code`      |
| `Name`         | Recommended | String |

### Entrez Email Requirement

For fetching sequences from NCBI databases (DNA/RNA), you need to provide a regular email address **or** an Entrez email for API access.

You can sign up for an account at [NCBI](https://www.ncbi.nlm.nih.gov/account/) and use your registered email address when running the program.

## af3map-chains
Creates a clear mapping between AlphaFold3 output chains and original input sequences for downstream analysis. Use after you converted the FASTA file to json with [AlphaFold3_tools](https://github.com/snufoodbiochem/Alphafold3_tools): fasta2json

### Usage
```bash
af3map-chains <input.fasta> <alphafold-output.json> [-o OUTPUT.tsv]
```
#### Example
```bash
af3map-chains af3_input.fasta prediction.json -o chain_map.tsv
```
---

### Example Workflow
1. Generate Input FASTA

```bash
af3build --email your@lab.com -o af3_input.fasta input_data.tsv
```

2. Create AlphaFold3 json
(using [AlphaFold3_tools | fasta2json](https://github.com/snufoodbiochem/Alphafold3_tools/blob/main/fasta2json.py))
```bash
python fasta2json.py af3_input.fasta
```

3. Generate Chain Mapping
```bash
af3map-chains -o chain_mapping.tsv af3_input.fasta af3_input.json
```

#### Example Output
| ChainID | ChainName           | CopyNumber | TotalCopies | SequenceType |
|---------|---------------------|------------|-------------|--------------|
| A       | Covid_SpikeProtein  | 1          | 2           | protein      |
| B       | Covid_SpikeProtein  | 2          | 2           | protein      |
| C       | Magnesium           | 1          | 1           | ligand       |
| D       | NM_001301244        | 1          | 3           | rna          |
| E       | NM_001301244        | 2          | 3           | rna          |
| F       | NM_001301244        | 3          | 3           | rna          |
| G       | CustomGeneName      | 1          | 2           | dna          |
| H       | CustomGeneName      | 2          | 2           | dna          |

---

### Validation Features
The af3map-chains command performs the following checks:

1. Chain Count Consistency
Ensures that the total number of chains in the JSON matches the expected number from the FASTA file (based on #N suffixes).

2. Order Preservation
Chains are mapped in the same order as they appear in the FASTA file.

3. Mismatch Handling
Reports mismatches between chain counts in the FASTA and JSON files and ensures no chains are missing or duplicated.


### PyMOL Integration
The generated mapping should be compatible to be used for downstream visualization and analysis in PyMOL eg.:

```python
import pandas as pd

# Load chain mapping
df = pd.read_csv("chain_mapping.tsv", sep="\t")

# Select all copies of SpikeProtein
spike_chains = ",".join(df[df["ChainName"] == "Covid_SpikeProtein"]["ChainID"])
cmd.select("spike", f"chain {spike_chains}")

# Color RNA chains
rna_chains = ",".join(df[df["SequenceType"] == "rna"]["ChainID"])
cmd.color("blue", f"chain {rna_chains}")
```

## af3estimate-tokens
You can estimate the number of tokens needed to run any Fasta file created by fasta2json.py. Based on this you can get a recommendation which GPU you should at least get to run your job.

> [!NOTE]
> You can use ```--smile_leniency``` to add a 5% headroom to your token count, as currently this tool cannot estimate the number of tokens from input other than dna,rna or protein

> [!NOTE]
> If you have tested or found more information on token numbers which successfully run on a different GPU please feel free to contact me or open a pull request to add to ```af3builder/gpu_config.py```

### Usage:
```bash
af3estimate-tokens [OPTIONS] <input.fasta>
```

#### Command-Line Options

| Option                | Description                                                                                   |
|-----------------------|-----------------------------------------------------------------------------------------------|
| `--verbose`           | Enable verbose output. Shows detailed information and breakdown for each sequence.            |
| `--recommendedgpu`    | Show GPU recommendations based on your input’s token count and current AlphaFold3 guidelines. |
| `--smile_leniency`    | Add 5% token headroom to the total token count to account for larger ligands/small molecules.<br>This helps ensure GPU recommendations remain robust for complex or large molecule inputs, as currently it does not calculate the number of tokens accurately from non DNA/RNA/Protein Sequences. |                                                            |


# Related Tools

- [AlphaFold3_tools](https://github.com/snufoodbiochem/Alphafold3_tools): Includes `fasta2json.py` for converting FASTA to JSON.
- [AlphaFold3 Official](https://github.com/google-deepmind/alphafold3): Main repository for AlphaFold3.

# Acknowledgements
This package builds around the work of [AlphaFold3_tools](https://github.com/snufoodbiochem/Alphafold3_tools) by Park et al.[^1] The methods and file formats provided here are adapted from and compatible with their toolkit. Please cite their publication if you use this package in your research:

This package also references code from [alphafold3](https://github.com/google-deepmind/alphafold3).

This package is intended to extend the current tools and provide a method for fast and high throughput creation of input for local installations of [alphafold3](https://github.com/google-deepmind/alphafold3).

[^1]: Park, Minho; Lee, Ungyu; Ha, Nam-Chul. "A rapid procedure to prepare for input file in a FASTA format for Alphafold3 local version." *Biodesign* 13, no. 1 (2025): 7–11.
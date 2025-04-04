import pandas as pd
import yaml
from pathlib import Path
from .exceptions import AF3Error
from .fetchers import SequenceFetcher
from .validator import ModificationValidator

class AF3Builder:
    def __init__(self, email=None):
        self.fetcher = SequenceFetcher(email)
        self.validator = ModificationValidator()

    def build(self, input_file, output_file):
        """Convert input file to AlphaFold3 FASTA"""
        df = self._read_input(input_file)
        fasta_lines = []

        if verbose == True:
            print("Input Dataframe:")
            print(df)

        for _, row in df.iterrows():
            entry = self._process_row(row)
            fasta_lines.append(entry)

        # Explicit empty lines between FASTA entries
        Path(output_file).write_text("\n\n".join(fasta_lines))

    def _read_input(self, input_file):
        """Read TSV/YAML input with strict validation"""
        path = Path(input_file)

        if path.suffix == ".tsv":
            # Read TSV with whitespace trimming and case normalization
            df = pd.read_csv(path, sep='\t')
            df.columns = df.columns.str.strip().str.upper()

            # Validate required columns
            required = {'ID', 'TYPE'}
            missing = required - set(df.columns)
            if missing:
                raise ValueError(
                    f"Missing columns in TSV: {missing}\n"
                    f"Found columns: {list(df.columns)}"
                )

            # Fix empty rows for Modifications
            df['MODIFICATIONS'] = df['MODIFICATIONS'].fillna('')

            return df.rename(columns={
                'TYPE': 'Type',
                'COPIES': 'Copies',
                'MODIFICATIONS': 'Modifications',
                'NAME': 'Name'
            })

        else:
            # Read and validate YAML
            with open(path) as f:
                data = yaml.safe_load(f)

            if not isinstance(data, list):
                raise ValueError("YAML must contain a list of entries")

            df = pd.DataFrame(data)

            # Normalize column names to UPPERCASE first
            df.columns = df.columns.str.strip().str.upper()

            # Check required columns exist
            required = {'ID', 'TYPE'}
            missing = required - set(df.columns)
            if missing:
                raise ValueError(
                    f"Missing fields in YAML: {missing}\n"
                    f"Found fields: {list(df.columns)}"
                )

            # Check for empty values
            if df[['ID', 'TYPE']].isnull().any().any():
                raise ValueError("YAML entries must have both 'ID' and 'TYPE'")

            # Fix empty rows for Modifications
            df['MODIFICATIONS'] = df['MODIFICATIONS'].fillna('')

            # Rename columns to match TSV format
            return df.rename(columns={
                'TYPE': 'Type',
                'COPIES': 'Copies',
                'MODIFICATIONS': 'Modifications',
                'NAME': 'Name'
            })


    def _process_row(self, row):
        """Handle one input row"""
        try:
            sequence = self.fetcher.get_sequence(row["ID"], row["Type"])

            if row["Type"].lower() == "rna":
                sequence = sequence.replace("T", "U")

            self.validator.validate(row.get("Modifications", ""), row["Type"])

            header = self._build_header(row, sequence_from_db=">" in sequence)
            return f"{header}\n{self._wrap_sequence(sequence)}"

        except AF3Error as e:
            return f"# ERROR: {str(e)}"

    def _build_header(self, row, sequence_from_db=False):
        """Ensure correct header order: name > mods > copies"""
        parts = []

        # Add database ID information if sequence is fetched from a database
        if sequence_from_db:
            parts.append(f">{row['ID']}")  # Original ID line from DB

        # Add custom name or molecule type prefix for FASTA2JSON compatibility
        seq_type = row['Type'].lower()
        name = row.get('Name', row['ID'])
        if seq_type in ['dna', 'rna', 'protein']:
            parts.append(f">{seq_type}{name}")
        elif seq_type in ['ligand', 'smile']:
            parts.append(f">{seq_type}{name}")

        if mods := row.get("Modifications"):
            parts[-1] += f" {mods}"  # Append modifications for FASTA2JSON to header

        if copies := row.get("Copies", 1) != 1:
            parts[-1] += f" #{copies}"  # Append oligomer count

        return "\n".join(parts)

    def _wrap_sequence(self, sequence):
        """Split long sequences into 60-character lines"""
        return "\n".join(sequence[i:i+60] for i in range(0, len(sequence), 60))

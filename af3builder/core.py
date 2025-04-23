import pandas as pd
import yaml
import re
from pathlib import Path
from .exceptions import AF3Error
from .fetchers import SequenceFetcher
from .validator import ModificationValidator

class AF3Builder:
    def __init__(self, email=None):
        self.fetcher = SequenceFetcher(email)
        self.validator = ModificationValidator()

    def build(self, input_file, output_file, verbose):
        """Convert input file to AlphaFold3 FASTA"""
        df = self._read_input(input_file)
        fasta_lines = []

        # Debugging
        if verbose:
            print("Input Dataframe:")
            print(df)

        # Split df into index and row data then process
        for _, row in df.iterrows():
            entry = self._process_row(row)
            fasta_lines.append(entry)

            # Debugging
            if verbose:
                print("Input Dataframe iteration:")
                print("row:")
                print(row)
                print("entry:")
                print(entry)
                print("fasta_lines:")
                print(fasta_lines)

        # Explicit empty lines between FASTA entries
        Path(output_file).write_text("\n\n".join(fasta_lines))

    def _read_input(self, input_file):
        """Read TSV/YAML input with accurate line-by-line preprocessing"""
        path = Path(input_file)

        if path.suffix == ".tsv":
            # TSV Handling
            df = pd.read_csv(
                path,
                sep='\t',
                dtype={'MODIFICATIONS': str, 'NAME': str},
                na_values=[''],
                keep_default_na=False
            )
            df.columns = df.columns.str.strip().str.upper()

            required = {'ID', 'TYPE'}
            if missing := required - set(df.columns):
                raise ValueError(f"Missing TSV columns: {missing}")

            df['MODIFICATIONS'] = df['MODIFICATIONS'].fillna('')
            df['NAME'] = df['NAME'].fillna('')

            return df.rename(columns={
                'TYPE': 'Type',
                'COPIES': 'Copies',
                'MODIFICATIONS': 'Modifications',
                'NAME': 'Name'
            })

        else:
            # YAML Line-by-Line Processing
            with open(path, 'r', encoding='utf-8') as f:
                yaml_lines = f.readlines()

            processed_lines = []
            for line in yaml_lines:
                # regex pattern for proper flag handling
                match = re.match(
                    r'^(\s*)(modifications|name)\s*:\s*(.*?)(?=\s*(#|$))',  # Removed inline flags
                    line,
                    flags=re.IGNORECASE  # Added flag parameter here
                )

                if match:
                    indent = match.group(1)
                    key = match.group(2).lower()
                    value = match.group(3).strip()
                    comment = re.search(r'\s*(#.*)$', line)
                    comment = comment.group(0) if comment else ''

                    # Handle empty values and special characters
                    needs_quotes = any([
                        not value,
                        any(c in value for c in {'&', '*', '!', '@'}),
                        ' ' in value,
                        value.startswith(('"', "'")),
                        value != value.strip()
                    ])

                    if needs_quotes:
                        processed_value = f'"{value}"' if value else '""'
                    else:
                        processed_value = value

                    new_line = f"{indent}{key}: {processed_value}{comment}\n"
                    processed_lines.append(new_line)
                else:
                    processed_lines.append(line)

            # Parse and validate YAML
            processed_yaml = ''.join(processed_lines)
            data = yaml.safe_load(processed_yaml)

            if not isinstance(data, list):
                raise ValueError("YAML must contain a list of entries")

            # Convert to DataFrame
            df = pd.DataFrame(data)
            df.columns = df.columns.str.strip().str.upper()

            # Validate required fields
            required = {'ID', 'TYPE'}
            if missing := required - set(df.columns):
                raise ValueError(f"Missing YAML fields: {missing}")

            # Ensure string types and handle missing values
            for col in {'MODIFICATIONS', 'NAME'} & set(df.columns):
                df[col] = df[col].fillna('').astype(str)

            return df.rename(columns={
                'TYPE': 'Type',
                'COPIES': 'Copies',
                'MODIFICATIONS': 'Modifications',
                'NAME': 'Name'
            })

    def _process_row(self, row):
        """Handle row with error catching"""
        try:
            original_header, sequence = self.fetcher.get_sequence(row["ID"], row["Type"])

            if row["Type"].lower() == "rna":
                sequence = sequence.replace("T", "U")

            self.validator.validate(row.get("Modifications", ""), row["Type"])
            return f"{self._build_header(row, original_header)}\n{self._wrap_sequence(sequence)}"

        except AF3Error as e:
            return f"# ERROR: {str(e)}"

    def _build_header(self, row, original_header=None):
        """Build FASTA header compatible with fasta2json's parser"""
        # Type Validation
        seq_type = row['Type'].lower()
        valid_types = {'protein', 'dna', 'rna', 'ligand', 'smile'}
        if seq_type not in valid_types:
            raise AF3Error(f"Invalid Type: {seq_type}. Allowed: {valid_types}")

        # Construct components in order: Name#Copies Modifications OriginalHeader
        components = []

        # Construct components
        name = row.get('Name', row['ID'])  # Use ID if Name is missing
        components = [name] if name else []

        # Modifications
        mods = row.get("Modifications", "").strip()
        if mods:
            components.append(mods)

        # Original header (strip existing > if present) (add | seperator)
        if original_header:
            components.append("|") # Seperator between DB FASTA Header and fasta2json modifications
            components.append(original_header.lstrip('>'))

        # Copies (always added as the last element)
        copies = int(row.get('Copies', 1))  # Default to 1 if missing or invalid
        components.append(f"#{copies}" if copies != 1 else "")  # Add only if >1


        return ">" + " ".join(components)

    def _wrap_sequence(self, sequence):
        """Split long sequences into 60-character lines"""
        return "\n".join(sequence[i:i+60] for i in range(0, len(sequence), 60))

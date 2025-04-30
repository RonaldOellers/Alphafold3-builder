import pandas as pd
import yaml
import re
from pathlib import Path
from .exceptions import AF3Error
from .fetchers import SequenceFetcher
from .validator import ModificationValidator
from .gpu_config import OFFICIAL_GPU_CAPACITY, COMMUNITY_GPU_CAPACITY

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
                dtype={'ID': str, 'TYPE':str, 'COPIES': int, 'MODIFICATIONS': str, 'NAME': str},
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
            
            # Ensure correct col types
            df['ID'] = df['ID'].astype(str)
            df['TYPE'] = df['TYPE'].astype(str)
            df['COPIES'] = df['COPIES'].astype(int)

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


###### Token and GPU Estimator ######
class AlphaFold3TokenCounter:
    """FASTA token counter with GPU recommendations based on AF3 docs"""
    
    def __init__(self, fasta_path, verbose=False, recommendedGPU=False, smile_leniency=False):
        self.fasta_path = fasta_path
        self.verbose = verbose
        self.recommendedGPU = recommendedGPU
        self.smile_leniency = smile_leniency
        self._token_data = {}
        self._original_tokens = 0

    def parse(self):
        current_header = None
        current_sequence = []
        current_multiplier = 1

        with open(self.fasta_path) as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Finalize previous sequence
                    if current_header:
                        seq = ''.join(current_sequence)
                        self._token_data[current_header] = len(seq) * current_multiplier
                    # Process new header
                    header_part = line[1:].split('#')
                    current_header = header_part[0].strip()
                    current_multiplier = int(header_part[1]) if len(header_part) > 1 and header_part[1].isdigit() else 1
                    current_sequence = []
                elif line:
                    current_sequence.append(line)
            # Add final sequence
            if current_header:
                seq = ''.join(current_sequence)
                self._token_data[current_header] = len(seq) * current_multiplier
        self._original_tokens = sum(self._token_data.values())

    @property
    def total_tokens(self):
        if self.smile_leniency: # adjust if user want headroom for smiles or ligands
            return int(self._original_tokens * 1.05)
        return self._original_tokens

    def _get_recommended_gpus(self):
        """Get GPUs that can handle token count"""
        adjusted_tokens = self.total_tokens
        valid_gpus = []
        
        # Check official GPUs
        for gpu, capacity in OFFICIAL_GPU_CAPACITY.items():
            if capacity >= adjusted_tokens:
                valid_gpus.append((gpu, capacity, False))
        
        # Check community GPUs
        for gpu, capacity in COMMUNITY_GPU_CAPACITY.items():
            if capacity >= adjusted_tokens:
                valid_gpus.append((gpu, capacity, True))
        
        return sorted(valid_gpus, key=lambda x: x[1])

    def summary(self):
        print(f"AlphaFold3 Token Report: {self.fasta_path}")
        print(f"Total Sequences: {len(self._token_data)}")
        
        if self.smile_leniency:
            print(f"Original Tokens: {self._original_tokens}")
            print(f"Adjusted Tokens (+5%): {self.total_tokens}\n")
        else:
            print(f"Total Tokens: {self.total_tokens}\n")

        if self.verbose:
            print("Sequence Breakdown:")
            for header, tokens in self._token_data.items():
                print(f"  {header}: {tokens} tokens")

        if self.recommendedGPU:
            print("\nHardware Recommendations:")
            valid_gpus = self._get_recommended_gpus()
            
            if valid_gpus:
                # Separate official and community GPUs
                official_gpus = [g for g in valid_gpus if not g[2]]
                community_gpus = [g for g in valid_gpus if g[2]]
                
                if official_gpus:
                    print("  Officially Supported:")
                    for gpu, cap, _ in official_gpus:
                        print(f"    - {gpu} ({cap} tokens)")
                        
                if community_gpus:
                    print("\n  Community-Reported:")
                    for gpu, cap, _ in community_gpus:
                        print(f"    - {gpu} ({cap} tokens)")
            else:
                print("  No supported GPUs can handle this input")
                print(f"  Largest supported: {max(OFFICIAL_GPU_CAPACITY.values())} tokens")
            
            print(40*'-')
            print("Please refer to https://github.com/google-deepmind/alphafold3/blob/main/docs/performance.md for official benchmarks and settings")
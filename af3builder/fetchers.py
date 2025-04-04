import requests
from Bio import Entrez
from .exceptions import SequenceFetchError, InvalidSequenceError

class SequenceFetcher:
    def __init__(self, email=None):
        self.email = email
        
    def get_sequence(self, identifier, seq_type):
        """Get sequence from ID or use raw input"""
        if self._is_raw_sequence(identifier):
            return self._clean_sequence(identifier, seq_type)
            
        seq_type = seq_type.lower()
        if seq_type == "protein":
            return self._fetch_uniprot(identifier)
        return self._fetch_ncbi(identifier)

    def _is_raw_sequence(self, identifier):
        """Check if input is a raw sequence (DNA, RNA, or protein) or a special case like ligand/smile."""
        # Handle special cases for ligands and SMILES strings
        if identifier.startswith(('ligand', 'smile')):
            return True

        # Define valid characters for DNA, RNA, and protein sequences
        valid_dna_rna_chars = set("ATUCG")  # DNA and RNA bases
        valid_protein_chars = set("ACDEFGHIKLMNPQRSTVWYUO")  # Amino Acids
        # Check if all characters are valid for DNA/RNA or protein
        identifier_upper = identifier.upper()
        if all(c in valid_dna_rna_chars for c in identifier_upper):
            return True
        if all(c in valid_protein_chars for c in identifier_upper):
            return True
        # If the identifier contains invalid characters (e.g., numbers, underscores), it's not a raw sequence
        return False
        
    def _clean_sequence(self, sequence, seq_type):
        """Remove spaces and validate"""
        clean_seq = sequence.upper().replace(" ", "")
        valid_chars = {
            "dna": "ATCGN",
            "rna": "AUCGN",
            "protein": "ACDEFGHIKLMNPQRSTVWY"
        }
        if not all(c in valid_chars[seq_type] for c in clean_seq):
            raise InvalidSequenceError(f"Invalid {seq_type} sequence")
        return clean_seq

    def _fetch_uniprot(self, uniprot_id):
        """Get protein sequence from UniProt"""
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
        try:
            response = requests.get(url)
            response.raise_for_status()
            return "".join(response.text.split("\n")[1:]).replace(" ", "")
        except Exception as e:
            raise SequenceFetchError(f"UniProt error: {str(e)}")

    def _fetch_ncbi(self, accession):
        """Get DNA/RNA sequence from NCBI"""
        if not self.email:
            raise ValueError("NCBI requires your email")
            
        Entrez.email = self.email
        try:
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta")
            record = handle.read()
            return "".join(record.split("\n")[1:]).replace(" ", "")
        except Exception as e:
            raise SequenceFetchError(f"NCBI error: {str(e)}")

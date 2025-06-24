import requests
from Bio import Entrez
from .exceptions import SequenceFetchError, InvalidSequenceError

class SequenceFetcher:
    def __init__(self, email=None):
        self.email = email

    def get_sequence(self, identifier, seq_type):
        """Get sequence from ID or use raw input"""
        if self._is_raw_sequence(identifier, seq_type):
            return None, self._clean_sequence(identifier, seq_type)  # (None, sequence)

        seq_type = seq_type.lower()
        if seq_type == "protein":
            return self._fetch_uniprot(identifier)
        return self._fetch_ncbi(identifier)

    def _is_raw_sequence(self, identifier, seq_type):
        """Check if input is raw sequence or special case"""
        seq_type = seq_type.lower()

        if not identifier or identifier.strip() == "":
            return False
        if seq_type in ['ligand', 'smile']:
            return True
        if identifier.startswith(('ligand', 'smile')):
            return True

        valid_dna_chars = set("ATCGN")
        valid_rna_chars = set("AUCGN")
        valid_protein_chars = set("ACDEFGHIKLMNPQRSTVWY")

        identifier_upper = identifier.upper()
        if seq_type == "dna":
            return all(c in valid_dna_chars for c in identifier_upper)
        elif seq_type == "rna":
            return all(c in valid_rna_chars for c in identifier_upper)
        elif seq_type == "protein":
            return all(c in valid_protein_chars for c in identifier_upper)
        return False

    def _clean_sequence(self, sequence, seq_type):
        """Clean and validate sequences"""
        clean_seq = sequence.upper().replace(" ", "")
        seq_type = seq_type.lower()

        if seq_type in ['ligand', 'smile']:
            return clean_seq  # Skip validation

        valid_chars = {
            "dna": "ATCGN",
            "rna": "AUCGN",
            "protein": "ACDEFGHIKLMNPQRSTVWY"
        }

        if not all(c in valid_chars[seq_type] for c in clean_seq):
            raise InvalidSequenceError(f"Invalid {seq_type} sequence")

        return clean_seq

    # In fetchers.py
    def _fetch_uniprot(self, uniprot_id):
        """Get UniProt entry with original header"""
        try:
            response = requests.get(f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta")
            response.raise_for_status()
            header, *seq_lines = response.text.splitlines()
            return header.lstrip('>'), ''.join(seq_lines).replace(' ', '')

        except Exception as e:
            raise SequenceFetchError(f"UniProt error: {str(e)}")


    def _fetch_ncbi(self, accession):
        """Get DNA/RNA sequence + header from NCBI"""
        if not self.email:
            raise ValueError("NCBI requires your email")

        Entrez.email = self.email
        try:
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta")
            record = handle.read().splitlines()
            return record[0].strip(), "".join(line.strip() for line in record[1:])
        except Exception as e:
            raise SequenceFetchError(f"NCBI error: {str(e)}")

import pytest
from af3builder.fetchers import SequenceFetcher
from af3builder.exceptions import SequenceFetchError, InvalidSequenceError
import time

# Rate limiting for NCBI requests (max 3 requests per second)
NCBI_DELAY = 1

def test_is_raw_sequence():
    fetcher = SequenceFetcher()

    # Valid DNA sequence
    assert fetcher._is_raw_sequence("ATGC", "dna") == True

    # Valid RNA sequence
    assert fetcher._is_raw_sequence("AUGCU", "rna") == True

    # Valid protein sequence
    assert fetcher._is_raw_sequence("MAEGASTERDA", "protein") == True

    # Ligand case
    assert fetcher._is_raw_sequence("ligand#3", "ligand") == True

    # Invalid identifier with numbers and underscores
    assert fetcher._is_raw_sequence("NM_001301244", "rna") == False

    # Invalid sequence with special characters
    assert fetcher._is_raw_sequence("ATGC_123", "dna") == False

    # Invalid sequence with non-biological characters
    assert fetcher._is_raw_sequence("XYZ12345", "protein") == False

    # Valid protein sequence with lowercase input (case-insensitivity)
    assert fetcher._is_raw_sequence("maegasterda", "protein") == True

    # Valid DNA sequence with lowercase input (case-insensitivity)
    assert fetcher._is_raw_sequence("atgc", "dna") == True

    # Empty input should return False
    assert fetcher._is_raw_sequence("", "protein") == False

    # Ligand prefix but invalid format should still return True
    assert fetcher._is_raw_sequence("ligand_invalid", "ligand") == True

def test_clean_sequence():
    fetcher = SequenceFetcher()

    # Valid DNA sequence
    assert fetcher._clean_sequence("ATGC", "dna") == "ATGC"

    # Valid RNA sequence
    assert fetcher._clean_sequence("AUGCU", "rna") == "AUGCU"

    # Valid protein sequence
    assert fetcher._clean_sequence("MAEGASTERDA", "protein") == "MAEGASTERDA"

    # Invalid DNA sequence
    with pytest.raises(InvalidSequenceError):
        fetcher._clean_sequence("ATGC_123", "dna")

    # Invalid RNA sequence
    with pytest.raises(InvalidSequenceError):
        fetcher._clean_sequence("AUGCU_123", "rna")

    # Invalid protein sequence
    with pytest.raises(InvalidSequenceError):
        fetcher._clean_sequence("XYZ12345", "protein")

def test_fetch_uniprot():
    fetcher = SequenceFetcher()

    # Test with known UniProt entry
    header, sequence = fetcher._fetch_uniprot("P0DTD1")

    assert "P0DTD1" in header
    assert len(sequence) > 1000  # Actual SARS-CoV-2 spike protein

    # Test with invalid ID
    with pytest.raises(SequenceFetchError):
        fetcher._fetch_uniprot("INVALID_ID_UNIPROT_123")

def test_fetch_ncbi():
    # NCBI requires email
    fetcher = SequenceFetcher(email="test@example.com")
    time.sleep(NCBI_DELAY)  # Rate limiting

    # Test with known NCBI RNA entry
    header, sequence = fetcher._fetch_ncbi("NM_001301244")

    assert "NM_001301244" in header
    assert len(sequence) > 100  # Actual mRNA sequence

    # Test with invalid accession
    time.sleep(NCBI_DELAY)
    with pytest.raises(SequenceFetchError):
        fetcher._fetch_ncbi("INVALID_ACCESSION_NCBI_123")

    # Test missing email
    fetcher_no_email = SequenceFetcher()
    with pytest.raises(ValueError):
        fetcher_no_email._fetch_ncbi("NM_001301244")

def test_get_sequence():
    fetcher = SequenceFetcher(email="test@example.com")

    # Test raw sequences
    header, sequence = fetcher.get_sequence("ATGC", "dna")
    assert header is None
    assert sequence == "ATGC"

    header, sequence = fetcher.get_sequence("augcuaaucgg", "rna")
    assert header is None
    assert sequence == "AUGCUAAUCGG"

    header, sequence = fetcher.get_sequence("maegasterda", "protein")
    assert header is None
    assert sequence == "MAEGASTERDA"

    # Test database lookups
    header, sequence = fetcher.get_sequence("P0DTD1", "protein")
    assert "P0DTD1" in header
    assert len(sequence) == 7096
    assert sequence.startswith("MESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVFIKRSDARTAPHGHVMVELVAELEGIQYGRS")

    time.sleep(NCBI_DELAY)

    header, sequence = fetcher.get_sequence("NM_001301244", "rna")
    assert "NM_001301244" in header
    assert len(sequence) == 1217
    assert sequence.startswith("GCTCGCACTCCCGCTCCTCCGCCCGACCGCGCGCTCGCCCCGCCGCTCCTGCTGCAGCCCCAGGGCCCCTCGCCGCCGCCACCATGGACGCCATCAAGAAGAAGATGCAGATGCTGA")

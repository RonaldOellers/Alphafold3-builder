import pytest  # Use pytest for running tests
from af3builder.fetchers import SequenceFetcher

def test_is_raw_sequence():
    fetcher = SequenceFetcher()

    # Valid DNA sequence
    assert fetcher._is_raw_sequence("ATGC") == True

    # Valid RNA sequence
    assert fetcher._is_raw_sequence("AUGCU") == True

    # Valid protein sequence
    assert fetcher._is_raw_sequence("MAEGASTERDA") == True

    # Ligand case
    assert fetcher._is_raw_sequence("ligand#3") == True

    # Invalid identifier with numbers and underscores
    assert fetcher._is_raw_sequence("NM_001301244") == False

    # Invalid sequence with special characters
    assert fetcher._is_raw_sequence("ATGC_123") == False

    # Invalid sequence with non-biological characters
    assert fetcher._is_raw_sequence("XYZ12345") == False

    # Valid protein sequence with lowercase input (case-insensitivity)
    assert fetcher._is_raw_sequence("maegasterda") == True

    # Valid DNA sequence with lowercase input (case-insensitivity)
    assert fetcher._is_raw_sequence("atgc") == True

    # Empty input should return False
    assert fetcher._is_raw_sequence("") == False

    # Ligand prefix but invalid format should still return True
    assert fetcher._is_raw_sequence("ligand_invalid") == True


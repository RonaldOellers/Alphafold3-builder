import pytest
from af3builder.fetchers import SequenceFetcher
from af3builder.exceptions import SequenceFetchError, InvalidSequenceError


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


def test_fetch_uniprot(mocker):
    fetcher = SequenceFetcher()

    # Mock a valid UniProt response
    mock_response = mocker.patch('requests.get')
    mock_response.return_value.status_code = 200
    mock_response.return_value.text = ">sp|P12345|Example\nMAEGASTERDA\nTTTGVCCSAQW"

    header, sequence = fetcher._fetch_uniprot("P12345")

    assert header == "sp|P12345|Example"
    assert sequence == "MAEGASTERDATTTVCCSAQW"

    # Mock an invalid UniProt response
    mock_response.return_value.status_code = 404
    with pytest.raises(SequenceFetchError):
        fetcher._fetch_uniprot("INVALID_ID")


def test_fetch_ncbi(mocker):
    fetcher = SequenceFetcher(email="test@example.com")

    # Mock a valid NCBI response
    mock_entrez = mocker.patch('Bio.Entrez.efetch')

    mock_handle = mock_entrez.return_value.__enter__.return_value  # Mock context manager behavior
    mock_handle.read.return_value = ">NM_001301244 Example RNA\nAUGCUAAUCGG\nCCGGAAUUCGG"

    header, sequence = fetcher._fetch_ncbi("NM_001301244")

    assert header == ">NM_001301244 Example RNA"
    assert sequence == "AUGCUAAUCGGCCGGAAUUCGG"

    # Missing email should raise ValueError
    fetcher_no_email = SequenceFetcher()

    with pytest.raises(ValueError):
        fetcher_no_email._fetch_ncbi("NM_001301244")

    # Mock an invalid NCBI response (e.g., network error)
    mock_entrez.side_effect = Exception("NCBI error")

    with pytest.raises(SequenceFetchError):
        fetcher._fetch_ncbi("INVALID_ACCESSION")


def test_get_sequence(mocker):
    fetcher = SequenceFetcher(email="test@example.com")

    # Mock raw DNA sequence handling
    header, sequence = fetcher.get_sequence("ATGC", "dna")

    assert header is None  # No original header for raw sequences
    assert sequence == "ATGC"

    # Mock raw RNA sequence handling (case-insensitive)
    header, sequence = fetcher.get_sequence("augcuaaucgg", "rna")

    assert header is None  # No original header for raw sequences
    assert sequence == "AUGCUAAUCGG"

    # Mock raw protein sequence handling (case-insensitive)
    header, sequence = fetcher.get_sequence("maegasterda", "protein")

    assert header is None  # No original header for raw sequences
    assert sequence == "MAEGASTERDA"

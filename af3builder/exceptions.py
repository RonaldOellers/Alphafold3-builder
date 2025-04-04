class AF3Error(Exception):
    """Base error for all custom errors"""
    
class SequenceFetchError(AF3Error):
    """Failed to get sequence from database"""
    
class InvalidModificationError(AF3Error):
    def __init__(self, code, seq_type):
        allowed = ", ".join(ALLOWED_MODS[seq_type])
        super().__init__(
            f"Invalid {seq_type} code: {code}\n"
            f"Allowed modifications: {allowed}"
        )
    
class InvalidSequenceError(AF3Error):
    """Bad DNA/RNA/Protein characters"""

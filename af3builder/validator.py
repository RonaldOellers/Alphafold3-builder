import re
from .exceptions import InvalidModificationError

class ModificationValidator:
    MOD_PATTERN = re.compile(r"&(\d+)_([A-Z0-9]{3})")
    ALLOWED_MODS = {
    "protein": [
        "3HY", "P1L", "PTR", "NMM", "LYZ", "CSO", "CSD", "OCS",
        "MYK", "SEP", "TPO", "NEP", "HIP", "ALY", "MLY", "M3L",
        "MLZ", "2MR", "AGM", "MCS", "HYP", "AHB", "SNN", "SNC",
        "TRF", "KCR", "CIR", "YHA"
    ],
    "dna": [
        "5CM", "C34", "5HC", "6OG", "6MA", "1CC", "8OG", "5FC", "3DR"
    ],
    "rna": [
        "PSU", "5MC", "OMC", "4OC", "5MU", "OMU", "UR3", "A2M",
        "MA6", "6MZ", "2MG", "OMG", "7MG", "RSQ"
    ]
}
    
    
    def validate(self, mod_string, seq_type):
        """Check modification syntax like &199_CSO"""
        if not mod_string:
            return
        
        seq_type = seq_type.lower()
        if seq_type not in ['protein', 'dna', 'rna']:
            raise ValueError(f"Invalid Type: {seq_type}")

        for match in self.MOD_PATTERN.finditer(mod_string):
            position, code = match.groups()
            
            if not position.isdigit():
                raise InvalidModificationError(f"Bad position: {position}")
                
            if code not in self.ALLOWED_MODS[seq_type.lower()]:
                raise InvalidModificationError(f"Invalid {seq_type} code: {code}")

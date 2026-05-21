from dataclasses import dataclass, field

@dataclass(slots=True)
class VariantTranscript:
    chromosome: str
    position: int
    reference: str
    alt: str
    cdna_transcript: str
    
    g_dot: str = None    
    gene: str = None
    c_dot: str = None
    p_dot1: str = None
    p_dot3: str = None
    protein_transcript: str = None
    strand: int = None
    splicing: str = None
    
    def __str__(self):
        return f"{self.chromosome}-{self.position}-{self.reference}-{self.alt} {self.cdna_transcript}:{self.c_dot}"

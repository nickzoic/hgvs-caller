"""Extract HGVS variant descriptions from a reference and a variant nucleotide sequence.

>>> from hgvs_extractor import extract_variant_string
>>> extract_variant_string("c.", "CATGACA", "CATATAACA", offset=2073)
'c.2077delinsATA'
"""

from .extractor import (
    TooManyVariationsException,
    extract_variant_string,
    extract_variants,
    reverse_complement,
    translate,
)

__version__ = "0.1.0"

__all__ = [
    "extract_variant_string",
    "extract_variants",
    "TooManyVariationsException",
    "reverse_complement",
    "translate",
    "__version__",
]

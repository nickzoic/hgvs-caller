"""Call HGVS variant descriptions from a reference and a variant nucleotide sequence.

>>> from hgvs_caller import find_variant_string
>>> find_variant_string("c.", "CATGACA", "CATATAACA", offset=2073)
'c.2077delinsATA'
"""

from .caller import (
    TooManyVariationsException,
    find_variant_string,
    find_variants,
    reverse_complement,
    translate,
)

__version__ = "0.1.0"

__all__ = [
    "find_variant_string",
    "find_variants",
    "TooManyVariationsException",
    "reverse_complement",
    "translate",
    "__version__",
]

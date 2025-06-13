#!/bin/env python
# -*- coding: utf-8 -*-
"""sequences.py
This module provides functions to help working with sequence data.
"""

import random

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


ALPHABETS = {
    "dna": list("ACGT"),
    "rna": list("ACGU"),
    "prot": list("ACDEFGHIKLMNPQRSTVWY"),
}


def generate_random_sequence(
    seqprefix: str = "seq",
    description: str = "Randomly-generated sequence",
    length: int = 1000000,
    alphabet="dna",
    seed=0,
) -> SeqRecord:
    """Returns a SeqRecord with random sequence of given length.
    
    This SeqRecord object also contains an operations attribute.
    This is an empty list by default.
    """
    if seed != 0:  # Allow for reproducible random sequence
        random.seed(seed)
    seqid = f"{seqprefix}_00000"
    record = SeqRecord(
        Seq("".join([random.choice(ALPHABETS[alphabet]) for _ in range(length)])),
        id=seqid,
        name=seqid,
        description=description,
    )
    record.operations = []  # Will hold structural operations
    return record

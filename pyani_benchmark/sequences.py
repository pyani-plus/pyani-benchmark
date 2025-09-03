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
    "at75": list("AAATTTCG"),
}


class PoolRecord(SeqRecord):
    """Subclass of Bio.SeqRecord

    The only real modification is the addition of an operations
    attribute, which holds structural changes to the genome sequence.
    """

    def __init__(self, *args, **kwargs) -> None:
        super(PoolRecord, self).__init__(*args, **kwargs)
        self.operations: list = []


def generate_random_sequence(
    seqprefix: str = "seq",
    description: str = "Randomly-generated sequence",
    length: int = 1000000,
    alphabet="dna",
    seed=0,
) -> PoolRecord:
    """Returns a PoolRecord with random sequence of given length.

    This PoolRecord object also contains an operations attribute.
    This is an empty list by default.
    """
    if seed != 0:  # Allow for reproducible random sequence
        random.seed(seed)
    seqid = f"{seqprefix}_00000"
    record = PoolRecord(
        Seq("".join([random.choice(ALPHABETS[alphabet]) for _ in range(length)])),
        id=seqid,
        name=seqid,
        description=description,
    )
    return record

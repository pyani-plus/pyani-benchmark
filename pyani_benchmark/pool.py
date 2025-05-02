#!/bin/env python
# -*- coding: utf-8 -*-
"""pool.py
This module contains the Pool class, which is used to manage a pool of
simulated genome sequences.
"""

from copy import deepcopy
from pathlib import Path

import networkx as nx

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class Pool:
    """A class representing a pool of simulated genome sequences."""

    def __init__(self, genome: SeqRecord, maxsize: int = 100, mutrate: float = 0.01):
        """Initialises the Pool with a seed genome."""
        self.__initialise__(genome, maxsize, mutrate)  # Initialise pool properties

    def __initialise__(self, genome, maxsize: int, mutrate: float):
        """Initialises empty pool data."""
        self._log: list[str] = ["Initialising pool..."]
        self._graph: nx.DiGraph = nx.DiGraph()  # Records relationships between genomes
        self._seqidx = 0  # Count of sequences, used to generate unique IDs
        self._pool: list = []  # Collection of genome SeqRecords

        # Constrain pool size/substitution rate
        self._maxsize = maxsize  # Maximum number of sequences in the pool
        self._mutrate = (
            mutrate  # Base mutation/substitution rate per base between generations
        )
        self._log.append(
            f"Initialised pool with maxsize {self._maxsize} and substitution rate {self._mutrate}"
        )

        # Add seed genome to pool
        self.__add_genome(genome)
        self._log.append(
            f"Initialised pool with genome {genome.id} (length: {len(genome.seq)})"
        )
        self._log.append("Pool initialised.")

    def __add_genome(self, genome):
        """Add a genome to the pool"""
        self._graph.add_node(genome.id)
        self._pool.append(genome)

    def write_log(self, fpath: Path):
        """Write the pool log to a file."""
        with fpath.open("w") as ofh:
            ofh.write("\n".join(self._log))

    def write_pool_dir(self, dirpath: Path):
        """Write individual pool sequences to a directory."""
        dirpath.mkdir(parents=True, exist_ok=True)
        for genome in self._pool:
            with (dirpath / f"{genome.id}.fasta").open("w") as ofh:
                SeqIO.write(genome, ofh, "fasta")

    def write_graph(self, fpath: Path):
        """Write the pool graph to a file."""
        nx.write_gml(self._graph, fpath)

    @property
    def log(self) -> list[str]:
        """Return a copy of the pool log."""
        return self._log[:]

    @property
    def pool(self) -> list[SeqRecord]:
        """Return a copy of the pool."""
        return deepcopy(self._pool)

    @property
    def graph(self) -> nx.DiGraph:
        """Return a copy of the pool graph."""
        return self._graph.copy()

    @property
    def difference_matrix(self):
        """Return the difference matrix for the pool."""
        return self._diffs.copy()

    def __len__(self) -> int:
        """Return pool size"""
        return len(self._pool)

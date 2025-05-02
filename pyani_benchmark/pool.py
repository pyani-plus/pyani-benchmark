#!/bin/env python
# -*- coding: utf-8 -*-
"""pool.py
This module contains the Pool class, which is used to manage a pool of
simulated genome sequences.
"""

import random

from copy import deepcopy
from pathlib import Path

import networkx as nx

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pyani_benchmark.sequences import ALPHABETS


class Pool:
    """A class representing a pool of simulated genome sequences."""

    def __init__(
        self,
        genome: SeqRecord,
        maxsize: int = 100,
        mutrate: float = 0.01,
        seqprefix="seq",
        alphabet="dna",
    ):
        """Initialises the Pool with a seed genome."""
        self.__initialise__(
            genome, maxsize, mutrate, seqprefix, alphabet
        )  # Initialise pool properties

    def __initialise__(
        self, genome, maxsize: int, mutrate: float, seqprefix: str, alphabet: str
    ):
        """Initialises empty pool data."""
        self._log: list[str] = ["Initialising pool..."]
        self._graph: nx.DiGraph = nx.DiGraph()  # Records relationships between genomes
        self._seqidx = 0  # Count of sequences, used to generate unique IDs
        self._pool: list = []  # Collection of genome SeqRecords
        self._seqprefix = seqprefix
        self._alphabet = ALPHABETS[alphabet]

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

    def __add_genome(self, genome) -> None:
        """Add a genome to the pool"""
        self._graph.add_node(genome.id)
        self._pool.append(genome)

    def mutate_genome(self, genome) -> SeqRecord:
        """Returns a copy of the passed genome, with symbol substitutions."""
        self._seqidx += 1  # Increment global sequence index
        new_seq = str(genome.seq)
        new_id = f"{self._seqprefix}+{self._seqidx:05d}"
        self._log.append(f"Generating genome {new_id} from genome {genome.id}")

        # Identify sequence positions to substitute symbols
        indices = [
            random.randrange(len(new_seq))
            for _ in range(round(len(new_seq) * self._mutrate))
        ]
        self._log.append(f"\tMutating {len(indices)} bases in genome {new_id}")

        # Make single-site substitutions
        for idx in indices:
            choices = list(set(self._alphabet) - set(new_seq[idx]))
            new_seq = new_seq[:idx] + random.choice(choices) + new_seq[idx + 1 :]

        return SeqRecord(
            Seq(new_seq), id=new_id, name=new_id, description="Evolved genome"
        )

    def evolve(self):
        """Update the pool by a single generation.

        We do this by generating a new updated sequence for each genome in the pool,
        shuffing the set of genomes, and then reducing the pool down to the maximum
        pool size, if necessary.
        """
        # Generate mutated genomes
        for genome in self._pool[:]:
            new_genome = self.mutate_genome(genome)  # Mutate genome
            self._log.append(
                f"Adding genome {new_genome.id} (length: {len(new_genome.seq)}) to the pool."
            )
            self._pool.append(new_genome)  # Update pool
            self._graph.add_edge(genome.id, new_genome.id)  # Update graph

        # Cull the pool back if necessary
        random.shuffle(self._pool)
        self._pool = self.pool[: self._maxsize]

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

#!/bin/env python
# -*- coding: utf-8 -*-
"""pool.py
This module contains the Pool class, which is used to manage a pool of
simulated genome sequences.
"""

import random

from copy import deepcopy
from functools import cache
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pyani_benchmark.sequences import ALPHABETS
from pyani_benchmark.utils import my_cmap


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
        self._startid = genome.id  # LUCA id

    @cache
    def __cached_difference_matrix(self, genomes):
        """Calculate the difference matrix for a set of genomes."""
        diffs = np.zeros((len(genomes), len(genomes)))
        for idx, genome in enumerate(genomes):
            for jdx, other in enumerate(genomes):
                diffs[idx, jdx] = self.__compare_genomes(genome, other)
        return diffs

    @cache
    def __compare_genomes(self, genome1, genome2):
        """Compare two genomes."""
        return sum([1 for a, b in zip(genome1, genome2) if a == b]) / len(genome1)

    def calc_difference_matrix(self):
        """Calculates the difference matrix for the pool."""
        seqpool = tuple([str(genome.seq) for genome in self._pool])
        self._diffs = self.__cached_difference_matrix(seqpool)

    def mutate_genome(self, genome, pointmutation=True) -> SeqRecord:
        """Returns a copy of the passed genome, with symbol substitutions."""
        self._seqidx += 1  # Increment global sequence index
        new_seq = str(genome.seq)
        new_id = f"{self._seqprefix}+{self._seqidx:05d}"
        self._log.append(f"Generating genome {new_id} from genome {genome.id}")

        # Make single-site substitutions, if point mutations are desired
        if pointmutation:
            # Identify sequence positions to substitute symbols
            indices = [
                random.randrange(len(new_seq))
                for _ in range(round(len(new_seq) * self._mutrate))
            ]
            self._log.append(f"\tMutating {len(indices)} bases in genome {new_id}")

            for idx in indices:
                choices = list(set(self._alphabet) - set(new_seq[idx]))
                new_seq = new_seq[:idx] + random.choice(choices) + new_seq[idx + 1:]

        return SeqRecord(
            Seq(new_seq), id=new_id, name=new_id, description="Evolved genome"
        )

    def evolve(self, pointmutation=True):
        """Update the pool by a single generation.

        We do this by generating a new updated sequence for each genome in the pool,
        shuffing the set of genomes, and then reducing the pool down to the maximum
        pool size, if necessary.
        """
        # Generate mutated genomes
        for genome in self._pool[:]:
            # Mutate genome or obtain direct copy, depending on setting
            new_genome = self.mutate_genome(genome, pointmutation)
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

    def write_pool(self, fpath):
        """Write the pool to a file."""
        SeqIO.write(self._pool, fpath, "fasta")

    def write_pool_dir(self, dirpath: Path):
        """Write individual pool sequences to a directory."""
        dirpath.mkdir(parents=True, exist_ok=True)
        for genome in self._pool:
            with (dirpath / f"{genome.id}.fasta").open("w") as ofh:
                SeqIO.write(genome, ofh, "fasta")

    def write_graph(self, fpath: Path):
        """Write the pool graph to a file."""
        nx.write_gml(self._graph, fpath)

    def draw_graph(self, fpath):
        """Draw the pool graph."""
        pool_names = [genome.id for genome in self._pool]
        colour_map = [
            "purple" if node in pool_names else "green" for node in self._graph
        ]
        plt.figure(1, figsize=(24, 24))
        nx.draw(
            self._graph,
            pos=nx.bfs_layout(self._graph, self._startid),
            # pos=nx.forceatlas2_layout(self._graph),
            with_labels=True,
            node_color=colour_map,
        )
        plt.savefig(fpath)
        plt.clf()

    def draw_heatmap(self, fpath):
        """Draw a heatmap of the difference matrix."""
        self.calc_difference_matrix()
        dfm = pd.DataFrame(self._diffs, index=[genome.id for genome in self._pool])
        dfm.columns = dfm.index
        sns.set(font_scale=0.5)
        sns.clustermap(dfm, cmap=my_cmap, vmin=0.8, vmax=1, figsize=(24, 24))
        plt.savefig(fpath)

    def write_difference_matrix(self, fpath):
        """Write the difference matrix for the pool."""
        self.calc_difference_matrix()
        np.savetxt(fpath, self._diffs, delimiter=",")

    def write_difference_dataframe(self, fpath):
        """Write the difference matrix for the pool as a dataframe."""
        self.calc_difference_matrix()
        dfm = pd.DataFrame(self._diffs, index=[genome.id for genome in self._pool])
        dfm.columns = dfm.index
        dfm.to_csv(fpath)

    def write_long_difference(self, fpath):
        """Write the differences to a file in long form."""
        with fpath.open("w") as ofh:
            dfm = pd.DataFrame(self._diffs, index=[genome.id for genome in self._pool])
            dfm.columns = dfm.index
            ofh.write("genome1,genome2,identity\n")
            for idx, row in dfm.iterrows():
                for jdx, value in row.items():
                    ofh.write(f"{idx},{jdx},{value}\n")

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

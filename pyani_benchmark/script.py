# -*- coding: utf-8 -*-
"""script.py
This code defines the CLI/script for the pyani-benchmark benchmarking tool.
"""

from pathlib import Path

import typer

from rich.progress import track, Progress, SpinnerColumn, TextColumn

from pyani_benchmark.cli_args import (
    Alphabet,
    OPT_ARG_TYPE_ALPHABET,
    OPT_ARG_TYPE_GENERATIONS,
    OPT_ARG_TYPE_POOLSIZE,
    OPT_ARG_TYPE_SEED,
    OPT_ARG_TYPE_SEQLEN,
    OPT_ARG_TYPE_SEQPREFIX,
    OPT_ARG_TYPE_SUBRATE,
    OPT_ARG_TYPE_OUTDIR,
    OPT_ARG_TYPE_VERSION,
)
from pyani_benchmark.pool import Pool
from pyani_benchmark.sequences import generate_random_sequence

app = typer.Typer()


@app.command()
def main(
    generations: OPT_ARG_TYPE_GENERATIONS = 10,
    seed: OPT_ARG_TYPE_SEED = 0,
    seqlen: OPT_ARG_TYPE_SEQLEN = 1000000,
    seqprefix: OPT_ARG_TYPE_SEQPREFIX = "seq",
    alphabet: OPT_ARG_TYPE_ALPHABET = Alphabet.dna,
    poolsize: OPT_ARG_TYPE_POOLSIZE = 100,
    subrate: OPT_ARG_TYPE_SUBRATE = 0.01,
    outdir: OPT_ARG_TYPE_OUTDIR = Path("./outdir"),
    version: OPT_ARG_TYPE_VERSION = False,
) -> None:
    """Entry point for the pyani-benchmark CLI."""
    # For now, we just print a message indicating that the script is running.
    # In a real implementation, you would call the main function of your benchmarking tool here.
    print("Running pyani-benchmark script...")
    print("Generating random sequence record...")
    record = generate_random_sequence(
        seqprefix=seqprefix, length=seqlen, alphabet=alphabet
    )
    print(f"\t...random sequence {record.id} has length {len(record)}.")
    print(f"Creating sequence pool, seeding with {record.id}....")
    # Generate a pool of synthetic genomes, with a single seed sequence
    pool = Pool(
        record,
        maxsize=poolsize,
        mutrate=subrate,
        seqprefix=seqprefix,
        alphabet=alphabet,
    )
    print(f"\t...sequence pool has length {len(pool)}.")
    print(f"\t{pool.pool[0]}")
    print(f"Evolving pool for {generations} generations.")
    # Iterate for the required number of generations, making the structural/sequence
    # changes specified at the CLI
    for generation in track(range(generations), description="Evolving..."):
        pool.evolve()
        print(f"Generation: {generation}, pool size: {len(pool)}")

    # At this point, the pairwise identity of the synthetic sequences
    # is determined: all columns in the sequence alignment represent a
    # base/amino acid with common evolutionary origin.
    # Calculate the pairwise distances and write these, and the underlying
    # genome pool, to an output directory.
    print(f"Calculating pairwise differences between {len(pool)} synthetic genomes.")
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        transient=False,
    ) as progress:
        progress.add_task(description="Calculating distances...", total=None)
        pool.calc_difference_matrix()
        progress.add_task(description="Writing intial pool sequences...", total=None)
        pool.write_pool_dir(outdir / "inital_pool" / "initial_pool_sequences")
        pool.write_pool(outdir / "inital_pool" / "initial_pool_sequences.fasta")
        progress.add_task(description="Writing intial pool graphs...", total=None)
        pool.write_graph(outdir / "inital_pool" / "inital_pool_graph.gml")
        pool.draw_graph(outdir / "inital_pool" / "inital_pool_graph.pdf")
        progress.add_task(description="Writing intial pool pairwise distances...", total=None)
        pool.write_difference_matrix(outdir / "inital_pool" / "inital_pool_distances.mat")
        pool.write_difference_dataframe(outdir / "inital_pool" / "inital_pool_distances.csv")
        pool.write_long_difference(outdir / "inital_pool" / "inital_pool_distances_long.csv")
        pool.draw_heatmap(outdir / "initial_pool" / "initial_pool_distances.pdf")

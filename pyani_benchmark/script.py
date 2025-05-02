# -*- coding: utf-8 -*-
"""script.py
This code defines the CLI/script for the pyani-benchmark benchmarking tool.
"""

import typer

from rich.progress import track

from pyani_benchmark.cli_args import (
    Alphabet,
    OPT_ARG_TYPE_ALPHABET,
    OPT_ARG_TYPE_GENERATIONS,
    OPT_ARG_TYPE_POOLSIZE,
    OPT_ARG_TYPE_SEED,
    OPT_ARG_TYPE_SEQLEN,
    OPT_ARG_TYPE_SEQPREFIX,
    OPT_ARG_TYPE_SUBRATE,
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
    for generation in track(range(generations), description="Evolving..."):
        pool.evolve()
        print(f"Generation: {generation}, pool size: {len(pool)}")

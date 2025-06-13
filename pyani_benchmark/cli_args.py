#!/bin/env python
# -*- coding: utf-8 -*-
"""cli_args.py
Defines command-line arguments in typer for the entry point.
"""

from enum import Enum
from pathlib import Path
from typing_extensions import Annotated

import typer

from pyani_benchmark import __app_name__, __app_version__


class Alphabet(str, Enum):
    """Defines permissible alphabet names for the CLI option"""

    dna = "dna"
    rna = "rna"
    prot = "prot"


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(f"{__app_name__} v{__app_version__}")
        raise typer.Exit()


OPT_ARG_TYPE_ALPHABET = Annotated[
    Alphabet, typer.Option("--alphabet", help="Seed sequence alphabet")
]

OPT_ARG_TYPE_GENERATIONS = Annotated[
    int,
    typer.Option("--generations", help="Number of generations to evolve genome pool."),
]

OPT_ARG_TYPE_SUBRATE = Annotated[
    float, typer.Option("--subrate", help="Per-symbol substitution probability.")
]

OPT_ARG_TYPE_INVRATE = Annotated[
    float, typer.Option("--invrate", help="Per-genome random inversion probability.")
]

OPT_ARG_TYPE_RECRATE = Annotated[
    float, typer.Option("--recrate", help="Per-genome recombination probability.")
]

OPT_ARG_TYPE_SHUFRATE = Annotated[
    float, typer.Option("--shufrate", help="Per-genome shuffling probability.")
]

OPT_ARG_TYPE_POOLSIZE = Annotated[
    int, typer.Option("--poolsize", help="Maximum size of evolving genome pool.")
]

OPT_ARG_TYPE_SEED = Annotated[
    int, typer.Option("--seed", help="Seed value for PRNG (0 sets no PRNG seed value).")
]

OPT_ARG_TYPE_SEQLEN = Annotated[
    int, typer.Option("--seqlen", help="Seed genome length.")
]

OPT_ARG_TYPE_SEQPREFIX = Annotated[
    str, typer.Option("--seqprefix", help="Prefix string for sequences.")
]

OPT_ARG_TYPE_OUTDIR = Annotated[
    Path, typer.Option("--outdir", help="Path to output directory.")
]

OPT_ARG_TYPE_VERSION = Annotated[
    bool,
    typer.Option(
        "--version",
        help="Show the application version and exit.",
        callback=_version_callback,
    ),
]

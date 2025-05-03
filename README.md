# `benchmark`

`benchmark` is a helper package for [`pyani-plus`](https://github.com/pyani-plus/pyani-plus) that generates mock bacterial genome sequences with known pairwise sequence identities. It is intended for benchmarking the accuracy and precision of genome distance calculation/estimation tools, including `pyani-plus`.

## Installation

1. Clone the repository from [https://github.com/pyani-plus/benchmark](https://github.com/pyani-plus/benchmark).
2. Change to the `benchmark` directory.
3. Install using `pip install -U -e .

The benchmark tool will be available as `pyani-benchmark`

```bash
git clone git@github.com:pyani-plus/benchmark.git
cd benchmark
pip install -U -e .
pyani-benchmark --help
```

## Usage

Invoke the tool with `pyani-benchmark`. To see inline help, use `pyani-benchmark --help`.

```bash
% pyani-benchmark --help                [10:39:45]

 Usage: pyani-benchmark [OPTIONS]

 Entry point for the pyani-benchmark CLI.


╭─ Options ───────────────────────────────────────────────────────────────────────────────────────────╮
│ --generations               INTEGER         Number of generations to evolve genome pool.            │
│                                             [default: 10]                                           │
│ --seed                      INTEGER         Seed value for PRNG (0 sets no PRNG seed value).        │
│                                             [default: 0]                                            │
│ --seqlen                    INTEGER         Seed genome length. [default: 1000000]                  │
│ --seqprefix                 TEXT            Prefix string for sequences. [default: seq]             │
│ --alphabet                  [dna|rna|prot]  Seed sequence alphabet [default: dna]                   │
│ --poolsize                  INTEGER         Maximum size of evolving genome pool. [default: 100]    │
│ --subrate                   FLOAT           Per-symbol substitution probability. [default: 0.01]    │
│ --version                                   Show the application version and exit.                  │
│ --install-completion                        Install completion for the current shell.               │
│ --show-completion                           Show completion for the current shell, to copy it or    │
│                                             customize the installation.                             │
│ --help                                      Show this message and exit.                             │
╰─────────────────────────────────────────────────────────────────────────────────────────────────────╯
```

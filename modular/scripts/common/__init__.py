"""Shared helpers for the modular analysis scripts.

Each analysis module (ARG, MGE, plasmid, taxonomy) is a standalone entry point.
Anything that was duplicated across the legacy scripts in ``scripts/`` lives here
so the modules stay independent of one another.

    maf.py     LAST .maf streaming and parsing, strand correction
    blocks.py  collapsing overlapping hits to the best one per region
    coords.py  fragment offsets and ``_partN`` -> read id
    tables.py  list-valued CSV columns, schema-preserving CSV writes
"""

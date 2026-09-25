#!/usr/bin/env bash
# Edit the root before submission. All database files stay in this location.
KRAKEN_ROOT="/path/to/kraken_storage"
KRAKEN_DB_NAME="cfrna_k2_bacteria_archaea_viral_human_fungi_full"
# Optional Unix group for shared access; empty uses normal ownership.
KRAKEN_SHARED_GROUP=""

KRAKEN_LIBRARIES=(archaea viral fungi human bacteria)
KRAKEN_KMER_LENGTH=35
KRAKEN_MINIMIZER_LENGTH=31
KRAKEN_MINIMIZER_SPACES=7

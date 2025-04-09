#!/bin/bash
f=${1##*/}
f=${f%.stin}
RIETAN=${0%/*}/STRUCTURE_TIDY
export RIETAN
SAMPLE_DIR=${1%/*}
export RIETAN
cd "$SAMPLE_DIR"
$RIETAN/structure_tidy $f.stin > $f-tmp.sto

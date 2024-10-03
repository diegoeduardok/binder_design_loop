#!/bin/bash

# Make sure you have edited your ~/.bashrc file to include colabfold_batch in your PATH. Or you can provide a full path to the executable.

colabfold_batch --num-recycle {NUM_RECYCLE} --amber --num-relax 1 --templates --custom-template-path {TEMPLATE_PATH} {FASTA_PATH} {OUTPUT_PATH}
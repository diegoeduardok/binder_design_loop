#!/bin/bash

PYTHON="/home/diego/miniconda3/envs/proteinmpnn/bin/python" # Add your own path
PROTEINMPNN="../ProteinMPNN/" # Add your own path

folder_with_pdbs="{INPUT_PDB_PATH}"

output_dir="{OUTPUT_PATH}"
if [ ! -d $output_dir ]
then
    mkdir -p $output_dir
fi


path_for_parsed_chains=$output_dir"/parsed_pdbs.jsonl"
path_for_assigned_chains=$output_dir"/assigned_pdbs.jsonl"
path_for_fixed_positions=$output_dir"/fixed_pdbs.jsonl"
chains_to_design="{CHAINS_TO_DESIGN}"
#The first amino acid in the chain corresponds to 1 and not PDB residues index for now.
design_only_positions="{POSITIONS_TO_DESIGN}" #design only these residues; use flag --specify_non_fixed

$PYTHON $PROTEINMPNN"/helper_scripts/parse_multiple_chains.py" --input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains

$PYTHON $PROTEINMPNN"/helper_scripts/assign_fixed_chains.py" --input_path=$path_for_parsed_chains --output_path=$path_for_assigned_chains --chain_list "$chains_to_design"

$PYTHON $PROTEINMPNN"/helper_scripts/make_fixed_positions_dict.py" --input_path=$path_for_parsed_chains --output_path=$path_for_fixed_positions --chain_list "$chains_to_design" --position_list "$design_only_positions" --specify_non_fixed

$PYTHON $PROTEINMPNN"/protein_mpnn_run.py" \
        --jsonl_path $path_for_parsed_chains \
        --chain_id_jsonl $path_for_assigned_chains \
        --fixed_positions_jsonl $path_for_fixed_positions \
        --out_folder $output_dir \
        --num_seq_per_target {SEQS_PER_TARGET} \
        --sampling_temp "0.1" \
        --seed 37 \
        --batch_size 1

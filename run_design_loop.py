'''NOTE: IF YOU ARE CHANGING FIXED RESIDUES OR CHAINS YOU NEED TO EDIT THIS FILE (LINES 137-138).
'''

import os
import shutil
import subprocess
from glob import glob
from natsort import natsorted
import numpy as np
from loss import compute_loss_designs

def run_colabfold(template_script, pdb_template_path, fasta_path, output_path, num_recycle=5):
    with open(template_script, "r") as infile:
        script = infile.read()
    script = script.format(
            TEMPLATE_PATH=pdb_template_path,
            FASTA_PATH=fasta_path,
            OUTPUT_PATH=output_path,
            NUM_RECYCLE=num_recycle,
        )
    if not os.path.isdir(output_path):
        os.makedirs(output_path, exist_ok=True)
    outpath_script = os.path.join(output_path, "command.sh")
    with open(outpath_script, "w") as outfile:
        outfile.write(script)
        os.chmod(outpath_script, 0o755) 
    subprocess.run(outpath_script, check=True)

def run_proteinMPNN(template_script, pdb_path, chains_to_design, positions_to_design, output_path, seqs_per_target):
    with open(template_script, "r") as infile:
        script = infile.read()
    script = script.format(
                INPUT_PDB_PATH=pdb_path,
                CHAINS_TO_DESIGN=chains_to_design,
                OUTPUT_PATH=output_path,
                POSITIONS_TO_DESIGN=positions_to_design,
                SEQS_PER_TARGET=seqs_per_target,
            )
    if not os.path.isdir(output_path):
        os.makedirs(output_path, exist_ok=True)
    outpath_script = os.path.join(output_path, "command.sh")
    with open(outpath_script, "w") as outfile:
        outfile.write(script)
        os.chmod(outpath_script, 0o755) 
    subprocess.run(outpath_script, check=True)

def prepare_sequence_to_structure(input_path, output_path, round_num):
    '''Based on ProteinMPNN output, prepare files for structure prediction.
    '''
    # Read FASTA file
    fasta_file = glob(os.path.join(input_path, "seqs", "*.fa"))
    assert(len(fasta_file) == 1)
    fasta_file = fasta_file[0]
    with open(fasta_file, "r") as infile:
        seqs = []
        for line in infile.readlines():
            if line[0] == ">":
                continue
            else:
                seqs.append(line.strip())

    # We need to append the constant CD20 sequence here
    cd20_seq = "MRESKTLGAVQIMNGLFHIALGGLLMIPAGIYAPICVTVWYPLWGGIMYIISGSLLAATEKNSRKCLVKGKMIMNSLSLFAAISGMILSIMDILNIKISHFLKMESLNFIRAHTPYINIYNCEPANPSEKNSPSTQYCYSIQSLFLGILSVMLIFAFFQELVIAG"

    # Write new file
    output_fasta = "" 
    for i, seq in enumerate(seqs[1:]):
        output_fasta += ">" + f"Round_{round_num}_seq_{i+1}\n"
        if i == len(seqs[1:]) - 1:
            output_fasta += cd20_seq + ":" + cd20_seq + ":" + seq
        else:
            output_fasta += cd20_seq + ":" + cd20_seq + ":" + seq + "\n"

    if not os.path.isdir(output_path):
        os.makedirs(output_path, exist_ok=True)

    with open(os.path.join(output_path, f"designed_sequences_round_{round_num}.fa"), "w") as outfile:
        outfile.write(output_fasta)

def prepare_structure_to_sequence(colabfold_results, output_path, round_num):
    '''Based on colabfold output, prepare files for sequence generation.
    '''
    # Compute loss to select best structure
    loss, pdb_paths = compute_loss_designs(colabfold_results)
    input_pdb_path = pdb_paths[np.argmin(loss)] # Best structure
    # Make sure file is there
    assert(os.path.isfile(input_pdb_path))

    return np.min(loss), input_pdb_path

def run_loop(outpath="./", initial_round=0, total_rounds=10):
    best_loss = np.inf
    chain_of_best_pdbs = []
    losses = []
    best_pdb = None
    # error_fasta_files = [] # These are files that returned errors with oc
    for r in range(initial_round, total_rounds):
        round_path = os.path.join(outpath, f"round_{r}")
        print(f"Starting round {r}. Will save results in {round_path}.")
        # Find fasta files and run colabfold
        fasta_path = os.path.join(outpath, f"round_{r}", f"designed_sequences_round_{r}.fa")
        colabfold_outpath = os.path.join(outpath, f"round_{r}", "colabfold_results")

        run_colabfold( # Need to put this in a try-except block...
            template_script="colabfold_template.sh",
            pdb_template_path="colabfold_pdb_templates/",
            fasta_path=fasta_path,
            output_path=colabfold_outpath,
        )
            
        # Prepare files for round r + 1
        output_path = os.path.join(outpath, f"round_{r+1}")
        min_batch_loss, best_batch_pdb = prepare_structure_to_sequence(colabfold_outpath, output_path, r+1)

        # Decide whether to accept designs or continue with previous round result
        if min_batch_loss < best_loss:
            chain_of_best_pdbs.append(best_batch_pdb)
            losses.append(min_batch_loss)
            best_loss = min_batch_loss
            best_pdb = best_batch_pdb
        else:
            print(min_batch_loss, "is greater than previous best loss:", best_loss)
            print("Using best model from previous round.")

        # Run sequence prediction
        protein_mpnn_output = os.path.join(outpath, f"round_{r+1}")

        if not os.path.isdir(protein_mpnn_output):
            os.makedirs(protein_mpnn_output, exist_ok=True)
            
        save_best_pdb = os.path.join(protein_mpnn_output, "previous_round.pdb")
        shutil.copy(os.path.abspath(best_pdb), save_best_pdb)
        
        run_proteinMPNN(
            template_script="proteinMPNN_template.sh", 
            pdb_path=protein_mpnn_output, 
            chains_to_design="C", 
            positions_to_design="2 3 4 6 8 9 10 11 12 13 14 15 16 17 18 19 20 21 25 26 27 28 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 72 74 75 77 79 80", 
            output_path=protein_mpnn_output, 
            seqs_per_target=30
        )

        prepare_sequence_to_structure(protein_mpnn_output, protein_mpnn_output, r+1)
        
        print(f"Finished round {r}. After this round, best PDB is {best_pdb}.")

    with open(os.path.join(outpath, "chain_of_best_models.csv"), "w") as outfile:
        oufile.write("model, loss" + "\n")
        for p, l in zip(chain_of_best_pdbs, losses):
            oufile.write(p + ", " + l + "\n")

def main():
    run_loop(outpath="./loop_oct_2/", initial_round=0, total_rounds=10)

if __name__ == "__main__":
    main()
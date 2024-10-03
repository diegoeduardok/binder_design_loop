# Binder design loop

Code to perform binder design loop on CD20 (lymphoma marker). The main script, `run_design_loop.py`, uses colabfold (i.e., AF2 Multimer) and ProteinMPNN sequentially and iteratively to optimize a peptide binder to the target protein.

## Installation

You should follow the instructions from [localcolabfold](https://github.com/YoshitakaMo/localcolabfold) to install that package. You will need to include the colabfold path in your PATH variable in `~/.bashrc` or edit the `colabfold_template.sh` file to include the full path to `colabfold_batch`.

Then, follow the instructions in [ProteinMPNN](https://github.com/dauparas/ProteinMPNN) to create a conda environment that can run those codes and clone the GitHub repo in your machine. You will need to provide the path to the Python executable from that environment and the paths to the ProteinMPNN codes in `proteinMPNN_template.sh`.

You will also need to have a access to a Python environment with `numpy`, `mdtraj`, and `natsort` installed to run the script `run_design_loop.py`. This environment can be the same as the one that runs ProteinMPNN. 

## Running the code

To run the code, first open `run_design_loop.py`. Edit line 153 according to your desired output path and total number of runs. 

```python
def main():
    run_loop(outpath="./loop_oct_2/", initial_round=0, total_rounds=10)
``` 

Notice that you can restart from a round larger than 0 if you already ran some of them, but the loop will have no memory of the best models produced in past rounds.

You can control the chains and residue positions (indexed from 1, not 0!) that will be designed by editing lines 137-138:

```python
run_proteinMPNN(
            template_script="proteinMPNN_template.sh", 
            pdb_path=protein_mpnn_output, 
            chains_to_design="C", # Modify as needed
            positions_to_design="2 3 4 6 8 9 10 11 12 13 14 15 16 17 18 19 20 21 25 26 27 28 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 72 74 75 77 79 80", # Modify as needed
            output_path=protein_mpnn_output, 
            seqs_per_target=30
        )
```

To start a loop, create the folder where the output will be saved. Then, in that directory, create a folder `round_0`. This folder must contain a fasta file named `designed_sequences_round_0.fa` to start the loop. If you have a PDB file with the initial guess for the binder and want to start from there, follow these steps to generate the fasta file:

1. Create a copy of `proteinMPNN_template.sh` (e.g., `proteinMPNN_init.sh`)
1. Edit `proteinMPNN_init.sh`: change `"{INPUT_PDB_PATH}"` with the path to the PDB file and change `"{OUTPUT_PATH}"` with the path to `round_0`.
1. Run `proteinMPNN_init.sh`. Check that ProteinMPNN created the output (.fa file) in the directory `round_0/seqs`.
1. Run the following code in an interactive Python session:

```python
>>>from run_design_loop import prepare_sequence_to_structure # May need to change import depending on where you stored the code
>>>prepare_sequence_to_structure("{output_path}/round_0", "{output_path}/round_0", 0) # Change output path here
```

After the last step, the file `round_0/designed_sequences_round_0.fa` should be available.

Now, you can run the desing loop:

```bash
$ python run_desing_loop.py
```


## Output

The output is organized as follows:

```
output_directory
	|
	|
	round_{r}
		|
		|_colabfold_results # Prdicted structures, scores, log file, etc.
		|
		|_seqs # Sequences designed by ProteinMPNN
		|
		|_designed_sequences_round_{r}.fa # Sequences from ProteinMPNN formatted for colabfold
		|
		|_Other files from ProteinMPNN...

```

## Known issues

For some unknown cause, `colabfold_batch` may produce segmentation faults after some rounds.

import os
from glob import glob
from natsort import natsorted
import json
import itertools
import numpy as np
import mdtraj as md

def compute_iPAE(pdb_file, pae_json, selection_1=None, selection_2=None):
    # Load PAE
    with open(pae_json, "r") as file:
        data = json.load(file)
        PAE = np.asarray(data['predicted_aligned_error'])
        # max_PAE = float(data['max_predicted_aligned_error'])

    # Find residue indices by atom index
    traj = md.load(pdb_file)
    res_indices = np.zeros((traj.n_atoms)).astype(int)
    for atom in traj.top.atoms:
        res_indices[atom.index] = atom.residue.index

    # Select atoms
    if selection_1 is None:
        selection_1 = "chainid 2 and name CB" # Resids may change
    if selection_2 is None:
        selection_2 = "chainid 0 1 and name CB and residue 100 to 140"
    A = traj.top.select(selection_1)
    B = traj.top.select(selection_2)
    pairs = np.asarray(list(itertools.product(A, B)))

    # Find contacts
    contacts = np.where((md.compute_distances(traj, pairs).reshape(A.shape[0], B.shape[0]) < 0.35) == True)
    selected_pairs = np.asarray([[res_indices[A[a]], res_indices[B[b]]] for a, b in zip(contacts[0],  contacts[1])])
    if len(selected_pairs) == 0:
        return np.inf
    total_selected_pairs = np.vstack([selected_pairs, selected_pairs[:, ::-1]]) # PAE matrix is not symmetric, must consider entry [a, b] and [b, a]
    
    # Find iPAE
    iPAE = np.median(PAE[total_selected_pairs[:, 0], total_selected_pairs[:, 1]])
    
    return iPAE


def read_score_file(score_file):
    with open(score_file, "r") as file:
        data = json.load(file)
        pLDDT = np.mean(np.asarray(data['plddt']))
        pTM = float(data['ptm'])
        ipTM = float(data['iptm'])
    return pLDDT, pTM, ipTM

def find_pLDDT_pTM_ipTM(score_files):
    results = []
    for score_file in score_files:
        pLDDT, pTM, ipTM = read_score_file(score_file)
        results.append([pLDDT, pTM, ipTM])
    return np.asarray(results)

def radius_of_g(pdb_file):
    traj = md.load(pdb_file)
    rg = md.compute_rg(traj.atom_slice(traj.top.select("chainid 2")))
    return rg[0]

def score(pLDDT, ipTM, iPAE, rg):
    '''This is the score to maximize.
    '''
    rg_penalty = np.zeros_like(rg)
    threshold = 2.5
    rg_penalty[np.where(rg > threshold)] = rg[np.where(rg > threshold)] - threshold
    return np.log(pLDDT / 100) + np.log(ipTM) - 0.2*np.log(iPAE) - rg_penalty # Make it less sensible to iPAE changes

def compute_loss_designs(dir_colabfold, selection_1=None, selection_2=None):
    # log_file = os.path.join(dir_colabfold, "log.txt")
    pdb_files = natsorted(glob(os.path.join(dir_colabfold, "*_*_relaxed_rank_001_alphafold2_multimer_v3_model_*_seed_000.pdb")))
    pae_files = natsorted(glob(os.path.join(dir_colabfold, "*_*_predicted_aligned_error_v1.json")))
    score_files = natsorted(glob(os.path.join(dir_colabfold, "*_*_scores_rank_001_alphafold2_multimer_v3_model_*_seed_000.json")))
    pLDDT_pTM_ipTM = find_pLDDT_pTM_ipTM(score_files)
    iPAE = np.asarray([compute_iPAE(pdb_file, pae_json, selection_1, selection_2) for pdb_file, pae_json in zip(pdb_files, pae_files)])
    rg = np.asarray([radius_of_g(pdb_file) for pdb_file in pdb_files])
    pLDDT = pLDDT_pTM_ipTM[:, 0]
    ipTM = pLDDT_pTM_ipTM[:, 2]
    L = -score(pLDDT, ipTM, iPAE, rg)
    return L, pdb_files
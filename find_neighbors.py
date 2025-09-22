import MDAnalysis as mda
import re
import numpy as np
from scipy.spatial import cKDTree as KDTree
import os
########################################################FUNCTIONS FOR FINDING CLOSEST RESIDUES
def find_closest_residues(universe, target_residue, cutoff=4, n=3):
    """
    Find the first n closest residues to a target residue within a given cutoff distance.

    Parameters:
    - universe: MDAnalysis Universe object.
    - target_residue: AtomGroup representing the target residue.
    - cutoff: float, distance cutoff in Å for the initial 'around' selection (default 10.0).
    - n: int, number of closest residues to return (default 5).

    Returns:
    - List of tuples (residue, distance) for the n closest residues.
    """
    # Step 1: Select residues around the target residue within the cutoff distance
    nearby_residues = universe.select_atoms(f"protein and around {cutoff} group target", target=target_residue)
    
    # Exclude the target residue itself
    nearby_residues = nearby_residues.residues - target_residue.residues

    # Step 2: Calculate distances between the center of mass of the target and nearby residues
    distances = []
    target_com = target_residue.center_of_mass()

    for residue in nearby_residues:
        distance = np.linalg.norm(target_com - residue.atoms.center_of_mass())
        distances.append((residue, distance))

    # Sort by distance and select the first n neighbors
    closest_residues = sorted(distances, key=lambda x: x[1])[:n]

    return closest_residues


def process_file_for_closest(input_file, output_file, universe, cutoff=5, n=3):
    """Read an input file, find closest neighbors, and write output with added columns."""
    protein = universe.select_atoms("protein")
    first_resid=protein.resids[0]
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        frame=0
        for line in infile:
            if "FRAME" in line:
                outfile.write(line+"\n")
                kdtree,residues,coms=build_kdtree(universe,frame)
                frame=frame+1
                continue
            parts = line.strip().split()
            elem=line.split()
            if len(elem) < 2:
                outfile.write("\n")                
                continue
            _, resname, resid = parts[0], parts[1],parts[2]
        
            # Select the target residue
            if "PAL" in line:
                target_residue = universe.select_atoms(f"resname {resname} and resid {resid} and name C1")
            else:
                target_residue = universe.select_atoms(f"resname {resname} and resid {resid}")
            if not target_residue:
                outfile.write(line.strip() + "\tNo_neighbors\n")
                continue
            
            # Find closest neighbors
            closest_residues = find_closest_residues_kdtree(target_residue,kdtree,residues,coms,resid,first_resid,coms,cutoff, n)
            resids=[]
            for resid in closest_residues:
                resids.append(universe.select_atoms(f"resid {resid}").residues[0])

            # Prepare output: add each closest neighbor as resname-resid
            neighbor_info = " ".join(f"{res.resname}{res.resid}" for res in resids)

            # Write the output line
            outfile.write(line.strip() + "\t" + neighbor_info + "\n")

def build_kdtree(universe, frame):
    """Build a KDTree for protein residues at a specific frame using their center of mass."""
    universe.trajectory[frame]
    
    # Select all protein atoms (not just CA atoms)
    protein = universe.select_atoms("protein")
    
    # Get the list of residues
    residues = protein.residues
    coms = []  # List to store the COMs of residues

    # Calculate the center of mass for each residue
    for residue in residues:
        atoms=residue.atoms
        com = atoms.center_of_mass()  # Center of mass of the residue
        coms.append(com)

    # Build the KDTree for the centers of mass of the residues
    kdtree = KDTree(coms)
    
    return kdtree, residues, coms

def find_closest_residues_kdtree(target_residue, kdtree, residues, positions,resid,first_resid,coms, cutoff=10.0, n=5):
    """Find the first n closest residues to a target residue using a precomputed KDTree."""
    target_com = target_residue.center_of_mass()
    indices = kdtree.query_ball_point(target_com, cutoff)
    indexes=[]
    distances=[]
    for idx in indices:
        distance = np.linalg.norm(target_com - coms[idx])
        distances.append(distance)
    if distances != []:
        closest_residues=[]
        paired=list(zip(indices,distances))
        sorted_pairs = sorted(paired, key=lambda x: x[1])
        indices = [x[0] for x in sorted_pairs]
        for idx in indices[:n]:
            closest_residues.append(idx+first_resid)
    else:
        closest_residues=[]

    return closest_residues

def find_closest_residues_kdtree_start_atom(target_residue, kdtree, residues, positions,resid,first_resid,coms, cutoff=10.0, n=5):
    """Find the first n closest residues to a target residue using a precomputed KDTree."""
    target_com = target_residue.positions[0]
    indices = kdtree.query_ball_point(target_com, cutoff)
    indexes=[]
    distances=[]
    for idx in indices:
        distance = np.linalg.norm(target_com - coms[idx])
        distances.append(distance)
    if distances != []:
        closest_residues=[]
        paired=list(zip(indices,distances))
        sorted_pairs = sorted(paired, key=lambda x: x[1])
        indices = [x[0] for x in sorted_pairs]
        for idx in indices[:n]:
            closest_residues.append(idx+first_resid)
    else:
        closest_residues=[]

    return closest_residues


########################################################FUNCTIONS FOR FINDING CLOSEST RESIDUES --END

def main():
    n=1 #number of closest neighbors
    cutoff=7 #cutoff to identify closest neighbors (same residue as within "cutoff" Angstroms)
    tolerance=0 #how many differences  can be in the channels in order to consider them to be the same
    topol="structure.gro"
    traj="traj.xtc"
    U=mda.Universe(topol,traj)
    filename="traj_path_nodes.dat"
    output="neighbors_"+filename
    if not os.path.exists(output):
        process_file_for_closest(filename, output, U, cutoff, n)

if __name__ == "__main__":
    main()


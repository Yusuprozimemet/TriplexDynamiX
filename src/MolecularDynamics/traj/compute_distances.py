import mdtraj as md
import numpy as np

def compute_distances(traj):
    # Print detailed topology information
    print("Detailed topology information:")
    for atom in traj.topology.atoms:
        print(f"Atom {atom.index}: {atom.name}, Residue {atom.residue.name}{atom.residue.index}, Chain {atom.residue.chain.index}")

    # Select atom indices for specific residues based on correct residue names
    residue_1_indices = traj.topology.select('resname C')
    residue_2_indices = traj.topology.select('resname A')
    residue_3_indices = traj.topology.select('resname U')

    # Print selected indices for debugging
    print(f"Residue 1 indices: {residue_1_indices}")
    print(f"Residue 2 indices: {residue_2_indices}")
    print(f"Residue 3 indices: {residue_3_indices}")

    # Check if indices arrays are empty
    if len(residue_1_indices) == 0:
        raise ValueError("No atoms found for residue C")
    if len(residue_2_indices) == 0:
        raise ValueError("No atoms found for residue A")
    if len(residue_3_indices) == 0:
        raise ValueError("No atoms found for residue U")

    # Compute center of mass for the selected residues
    cen1 = md.compute_center_of_mass(traj.atom_slice(residue_1_indices))
    cen2 = md.compute_center_of_mass(traj.atom_slice(residue_2_indices))
    cen3 = md.compute_center_of_mass(traj.atom_slice(residue_3_indices))

    # Calculate distances between centers of mass
    squared_dist1 = np.sum(np.square(cen1 - cen2), axis=1)
    dist1 = np.sqrt(squared_dist1)
    squared_dist2 = np.sum(np.square(cen1 - cen3), axis=1)
    dist2 = np.sqrt(squared_dist2)
    squared_dist3 = np.sum(np.square(cen2 - cen3), axis=1)
    dist3 = np.sqrt(squared_dist3)
    s = (dist1 + dist2 + dist3) / 2
    area = np.sqrt(s * (s - dist1) * (s - dist2) * (s - dist3))

    return dist1, dist2, dist3, area

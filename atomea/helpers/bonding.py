# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

import numpy as np
import numpy.typing as npt
from scipy.sparse import csr_array
from scipy.sparse.csgraph import connected_components


def get_molecule_ids(bonds: npt.NDArray[np.uint64], n_atoms: int):
    """Deterministic
    [molecule IDs][containers.ensemble.topology.ids.IDs.molecules] based on
    [covalent bonding][containers.ensemble.topology.connectivity.bonds.Bonds].
    IDs are assigned based on the following rules:

    - The largest molecules get the lowest IDs.
    - Tie breakers for atoms with the same number of atoms is done by the
        lowest atom index contained within the molecule getting lower IDs.

    Args:
        bonds: A 2D numpy array where each row represents a bond
            [atom_index_0, atom_index_1].

    Returns:
        A 1D numpy array where each element at index `i` is the
            deterministic id_molecule for atom `i`.
    """
    assert isinstance(bonds, np.ndarray)
    assert bonds.ndim == 2
    assert bonds.shape[1] == 2
    assert n_atoms > 0

    # Create an adjacency matrix
    row_indices = bonds[:, 0]
    col_indices = bonds[:, 1]

    data = np.ones(len(bonds) * 2)
    adj_array = csr_array(
        (
            data,
            (
                np.concatenate([row_indices, col_indices]),
                np.concatenate([col_indices, row_indices]),
            ),
        ),
        shape=(n_atoms, n_atoms),
    )

    # Find connected components (initial molecule IDs)
    n_components, initial_labels = connected_components(
        csgraph=adj_array, directed=False, return_labels=True
    )

    # Initialize final molecule IDs to -1 for all atoms
    final_id_molecule = np.full(n_atoms, -1, dtype=np.int64)

    # Gather properties for each distinct molecule found by connected_components
    molecule_properties = {}
    for initial_mol_id in range(n_components):
        atom_indices_in_mol = np.where(initial_labels == initial_mol_id)[0]

        if len(atom_indices_in_mol) > 0:
            mol_size = len(atom_indices_in_mol)
            min_atom_idx = np.min(atom_indices_in_mol)

            molecule_properties[initial_mol_id] = {
                "size": mol_size,
                "min_atom_idx": min_atom_idx,
                "atom_indices": atom_indices_in_mol,
            }

    # Handle isolated atoms (those not in any bond but within total_n_atoms scope)
    # These will have an initial_label of -1 from connected_components, or
    # If this atom index hasn't been assigned an ID (i.e., it's still -1 from initialization) AND
    next_isolated_mol_id = n_components  # Start isolated IDs after bonded molecules
    for atom_idx in range(n_atoms):
        # If this atom index hasn't been assigned an ID (i.e., it's still -1 from initialization) AND
        # it was NOT found in any bond (meaning it's truly isolated).
        if final_id_molecule[atom_idx] == -1:
            # Treat each isolated atom as a molecule of size 1
            isolated_atom_mol_id = next_isolated_mol_id
            molecule_properties[isolated_atom_mol_id] = {
                "size": 1,
                "min_atom_idx": atom_idx,  # For isolated atoms, their own index is the min
                "atom_indices": np.array([atom_idx]),
            }
            next_isolated_mol_id += 1

    # Sort molecules based on the specified rules:
    # 1. Largest size first (descending: -size)
    # 2. Then by lowest atom index (ascending: min_atom_idx)
    sortable_molecules = []
    for initial_mol_id, props in molecule_properties.items():
        size = props["size"]
        min_atom_idx = props["min_atom_idx"]
        sortable_molecules.append((initial_mol_id, size, int(min_atom_idx)))

    sortable_molecules.sort(key=lambda x: (-x[1], x[2]))

    # Assign new deterministic molecule_ids
    new_mol_id_counter = 0
    initial_to_new_id_map = {}

    for initial_mol_id, _, _ in sortable_molecules:
        initial_to_new_id_map[initial_mol_id] = new_mol_id_counter
        new_mol_id_counter += 1

    for initial_mol_id, props in molecule_properties.items():
        if initial_mol_id in initial_to_new_id_map:
            new_id = initial_to_new_id_map[initial_mol_id]
            idxs = np.where(initial_labels == initial_mol_id)[0]
            final_id_molecule[idxs] = new_id

    return final_id_molecule

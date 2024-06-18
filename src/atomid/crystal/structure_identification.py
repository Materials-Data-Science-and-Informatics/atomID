"""Crystal structure identification and lattice parameter calculation."""

from math import sqrt
from typing import Tuple

import numpy as np


def get_crsytal_structure_from_id(id: int) -> str:
    """
    Get the crystal structure type from the given ID.

    Parameters
    ----------
    id : int
        The ID of the crystal structure.

    Returns
    -------
    str
        The crystal structure type.
    """
    structure_type = {
        0: "other",
        1: "fcc",
        2: "hcp",
        3: "bcc",
        4: "ico",
        5: "sc",
        6: "cubic diamond",
        7: "hex diamond",
        8: "graphene",
    }
    return structure_type.get(id, "other")


def analyse_polyhedral_template_matching_data(
    atoms_structure_type: np.ndarray,
) -> Tuple:
    """Analyse polyhedral template data to identify crystal structure."""
    unique, counts = np.unique(atoms_structure_type, return_counts=True)
    structure_id = unique[np.argmax(counts)]
    return structure_id, get_crsytal_structure_from_id(structure_id)


def find_lattice_parameter(
    interatomic_distance: np.ndarray,
    structure_type_atoms: np.ndarray,
    crystal_type_id: int,
) -> float:
    # create mask for the structure type matching the crystal type id
    mask = structure_type_atoms == crystal_type_id

    # get the interatomic distances for the given structure type
    filtered_interatomic_distance = interatomic_distance[mask]

    # calculate mean of the interatomic distances
    mean_interactomic_distance = np.mean(filtered_interatomic_distance)

    # calculate multiplier for the lattice parameter

    multiplier = {
        1: sqrt(2),
        2: 1,
        3: sqrt(4 / 3),
        4: 4 / sqrt(3),
        5: 1,
        6: sqrt(16 / 3),
        7: sqrt(16 / 3),
        8: sqrt(2),
    }

    lattice_parameter: float = round(
        mean_interactomic_distance * multiplier[crystal_type_id], 3
    )
    return lattice_parameter

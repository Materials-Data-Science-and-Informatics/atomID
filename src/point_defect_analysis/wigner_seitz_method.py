"""Point defect identification using the Wigner-Seitz method."""

from typing import Tuple, Optional

import numpy as np


def find_nearest_atom(atom: tuple, atom_positions: np.ndarray) -> Tuple[int, list]:
    """Find the nearest atom to a given defect position.

    Parameters
    ----------
    atom : tuple
        The position of the defect atom.
    atom_positions : np.ndarra
        The positions of the atoms in the lattice.

    Returns
    -------
    nearest_index : int
        The index of the nearest atom.
    distance : float
        The distance between the defect and the nearest atom.
    """
    distances = np.linalg.norm(atom_positions - atom, axis=1)
    nearest_index: int = np.argmin(distances)
    return nearest_index, distances[nearest_index]


def analyze_defects(
    reference_positions: list,
    actual_positions: list,
    species_ref: Optional[list] = None,
    species_actual: Optional[list] = None,
) -> dict:
    """Analyze the lattice for vacancy and interstitial defects.

    Parameters
    ----------
    reference_positions : list of tuples
        The expected positions of the atoms in the lattice.
    actual_positions : list of tuples
        The actual positions of the atoms in the lattice.

    Returns
    -------
    defect_analysis : dict
        A dictionary containing the vacancy and interstitial defects.

    """
    reference_positions = np.array(reference_positions)
    actual_positions = np.array(actual_positions)
    atom_position_count = np.zeros(len(reference_positions))
    substitution_count = np.zeros(len(reference_positions))

    # Process actual positions and compare with reference to identify defects
    for i, actual in enumerate(actual_positions):
        nearest_index, _ = find_nearest_atom(actual, reference_positions)
        atom_position_count[nearest_index] += 1

        # Check for substitutions if species information is provided
        if species_actual and species_ref:
            if species_actual[i] != species_ref[nearest_index]:
                substitution_count[nearest_index] += 1

    # Determine vacancies taking into account both atom positions and substitutions
    vacancies = [
        (i, tuple(pos))
        for i, pos in enumerate(reference_positions)
        if atom_position_count[i] == 0
        and (not species_actual or substitution_count[i] == 0)
    ]

    vacancies = [
        (i, tuple(pos))
        for i, pos in enumerate(reference_positions)
        if atom_position_count[i] == 0
    ]
    interstitials = [
        (i, tuple(pos))
        for i, pos in enumerate(actual_positions)
        if atom_position_count[i] > 1
    ]
    substitutions = [
        (i, tuple(pos))
        for i, pos in enumerate(reference_positions)
        if substitution_count[i] > 0
    ]

    return {
        "Vacancies": {
            "count": len(vacancies),
            "fraction": len(vacancies) / len(reference_positions),
        },
        "Interstitials": {
            "count": len(interstitials),
            "fraction": len(interstitials) / len(actual_positions),
        },
        "Substitutions": {
            "count": len(substitutions),
            "fraction": len(substitutions) / len(reference_positions),
        },
    }

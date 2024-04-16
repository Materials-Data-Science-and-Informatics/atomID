import numpy as np
from typing import Tuple


def find_nearest_atom(atom: tuple, atom_positions: np.ndarray) -> Tuple[int, float]:
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
    nearest_index = np.argmin(distances)
    return nearest_index, distances[nearest_index]


def analyze_defects(
    reference_positions: list,
    actual_positions: list,
    species_ref: list = [],
    species_actual: list = [],
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
    atom_position_list = np.zeros(len(reference_positions))
    substitution_list = np.zeros(len(reference_positions))
    # Analyze atom positions to detect vacancies and interstitials

    if len(species_actual) != 0:
        for i, actual in enumerate(actual_positions):
            nearest_index, _ = find_nearest_atom(actual, reference_positions)
            if species_actual[nearest_index] == species_ref:
                if species_actual[i] == species_ref[nearest_index]:
                    atom_position_list[nearest_index] += 1
                else:
                    substitution_list[nearest_index] += 1
        vacancies = [
            (i, tuple(pos))
            for i, pos in enumerate(reference_positions)
            if atom_position_list[i] == 0 and substitution_list[i] == 0
        ]

    else:
        for actual in actual_positions:
            nearest_index, _ = find_nearest_atom(actual, reference_positions)
            atom_position_list[nearest_index] += 1
        vacancies = [
            (i, tuple(pos))
            for i, pos in enumerate(reference_positions)
            if atom_position_list[i] == 0
        ]

    vacancy_count = len(vacancies)
    vacancy_fraction = round(vacancy_count / len(actual_positions), 3)
    interstitials = [
        (i, tuple(pos))
        for i, pos in enumerate(actual_positions)
        if atom_position_list[i] > 1
    ]
    interstitial_count = len(interstitials)
    interstitial_fraction = round(interstitial_count / len(actual_positions), 3)

    substitutions = [
        (i, tuple(pos))
        for i, pos in enumerate(reference_positions)
        if substitution_list[i] > 0
    ]
    substitution_count = len(substitutions)
    substitution_fraction = round(substitution_count / len(actual_positions), 3)
    return {
        "Vacancies": {"count": vacancy_count, "fraction": vacancy_fraction},
        "Interstitials": {
            "count": interstitial_count,
            "fraction": interstitial_fraction,
        },
        "Substitutions": {
            "count": substitution_count,
            "fraction": substitution_fraction,
        },
    }

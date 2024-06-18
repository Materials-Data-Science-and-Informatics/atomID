"""Wigner-Seitz method for point defect analysis."""

from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
from annoy import AnnoyIndex
from sklearn.neighbors import KDTree


def analyze_defects(
    reference_positions: List[Tuple[float, float, float]],
    actual_positions: List[Tuple[float, float, float]],
    species_ref: Optional[List[str]] = None,
    species_actual: Optional[List[str]] = None,
    method: Optional[str] = None,
) -> Dict[str, Dict[str, float]]:
    """Analyze the lattice for vacancy, interstitial, and substitution defects.

    Parameters
    ----------
    reference_positions : list of tuples
        The expected positions of the atoms in the lattice.
    actual_positions : list of tuples
        The actual positions of the atoms in the lattice.
    species_ref : list of str, optional
        Species at each reference position.
    species_actual : list of str, optional
        Species at each actual position.
    method : str, optional
        The method to find nearest positions ('annoy' for using AnnoyIndex).

    Returns
    -------
    dict
        A dictionary containing the counts and concentration of vacancies, interstitials, and substitutions.
    """
    reference_array = np.array(reference_positions)
    actual_array = np.array(actual_positions)

    atom_position_count = np.zeros(len(reference_array))
    substitution_count = np.zeros(len(reference_array))
    index_finder = create_index_finder(reference_array, method)

    if species_ref is not None and species_actual is not None:
        identify_substitution = True
    else:
        identify_substitution = False

    for i, actual in enumerate(actual_array):
        nearest_index = index_finder(actual)
        atom_position_count[nearest_index] += 1
        if identify_substitution and (species_actual[i] != species_ref[nearest_index]):
            substitution_count[nearest_index] += 1

    defects: dict = calculate_defects(
        reference_array, atom_position_count, substitution_count
    )
    return defects


def create_index_finder(
    reference_array: np.ndarray, method: Optional[str] = None
) -> Callable:
    """Create a function to find the index of the nearest reference position.

    Parameters
    ----------
    reference_array : np.ndarray
        The reference positions of the atoms.
    method : str, optional
        The method to find nearest positions ('annoy' for using AnnoyIndex).

    Returns
    -------
    function
    A function that takes an actual position and returns the index of the nearest reference position.
    """
    if method == "annoy":
        t = AnnoyIndex(len(reference_array[0]), "euclidean")
        for i, ref in enumerate(reference_array):
            t.add_item(i, ref)
        t.build(10)
        return lambda x: t.get_nns_by_vector(x, 1)[0]
    elif method == "kd_tree":
        kdtree = KDTree(reference_array)
        return lambda x: kdtree.query(x.reshape(1, -1), k=1)[1][0][0]
    else:
        return lambda x: np.argmin(np.sum((reference_array - x) ** 2, axis=1))


def calculate_defects(
    reference_array: np.ndarray,
    atom_position_count: np.ndarray,
    substitution_count: np.ndarray,
) -> Dict:
    """Calculate the number and concentration of vacancies, interstitials, and substitutions.

    Parameters
    ----------
    reference_array : np.ndarray
        The reference positions of the atoms.
    atom_position_count : np.ndarray
        The number of atoms at each reference position.
    substitution_count : np.ndarray
        The number of substitutions at each reference position.

    Returns
    -------
    dict
        A dictionary containing the counts and concentration of vacancies, interstitials, and substitutions.
    """
    vacancies = [
        (i, tuple(pos))
        for i, pos in enumerate(reference_array)
        if atom_position_count[i] == 0
    ]
    interstitials = [
        (i, tuple(pos))
        for i, pos in enumerate(reference_array)
        if atom_position_count[i] > 1
    ]
    substitutions = [
        (i, tuple(pos))
        for i, pos in enumerate(reference_array)
        if substitution_count[i] > 0
    ]

    return {
        "vacancies": {
            "count": len(vacancies),
            "concentration": len(vacancies) / len(reference_array),
        },
        "interstitials": {
            "count": len(interstitials),
            "concentration": len(interstitials) / len(reference_array),
        },
        "substitutions": {
            "count": len(substitutions) - len(interstitials),
            "concentration": (len(substitutions) - len(interstitials))
            / len(reference_array),
        },
    }

"""Crystal structure identification and lattice parameter calculation."""

import logging
from math import sqrt
from typing import Tuple

import numpy as np
from atomrdf import System
from scipy.signal import find_peaks


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


def get_crystal_structure_using_cna(pyscal_system: System) -> str:
    """
    Get the crystal structure using adaptive common neighbour analysis.

    Parameters
    ----------
    pyscal_system : pyscal.System
        The pyscal system object.

    Returns
    -------
    str
        The identified crystal structure type.
    """
    cna_results = pyscal_system.analyze.common_neighbor_analysis()
    logging.info("Adaptive common neighbour analysis results: %s", cna_results)

    # Find the most frequent crystal structure from CNA results.
    crystal_type = max(cna_results, key=cna_results.get)

    if crystal_type == "others":
        # Further analyse if the crystal type could be diamond related.
        crystal_type = analyse_diamond_structures(pyscal_system)

    logging.info("Selected crystal structure type: %s", crystal_type)
    return str(crystal_type)


def analyse_diamond_structures(pyscal_system: System) -> str:
    """
    Analyse diamond structures and identify the dominant type if any.

    Parameters
    ----------
    pyscal_system : pyscal.System
        The pyscal system object.

    Returns
    -------
    str
        The dominant diamond structure type or 'others' if none found.
    """
    diamond_results = pyscal_system.analyze.diamond_structure()
    logging.info("Initial diamond structure analysis results: %s", diamond_results)

    # Aggregate the diamond structure counts.
    diamond_results["cubic diamond"] = sum(
        diamond_results.get(key, 0)
        for key in ["cubic diamond", "cubic diamond 1NN", "cubic diamond 2NN"]
    )
    diamond_results["hex diamond"] = sum(
        diamond_results.get(key, 0)
        for key in ["hex diamond", "hex diamond 1NN", "hex diamond 2NN"]
    )

    # Clean up the dictionary.
    for key in [
        "cubic diamond 1NN",
        "cubic diamond 2NN",
        "hex diamond 1NN",
        "hex diamond 2NN",
    ]:
        diamond_results.pop(key, None)

    logging.info("Consolidated diamond structure analysis results: %s", diamond_results)

    # Choose the most frequent diamond structure.
    diamond_type = max(diamond_results, key=diamond_results.get, default="others")
    return str(diamond_type)


def find_lattice_parameter_2(
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
        2: sqrt(2),
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


def find_lattice_parameter(crystal_system: System, lattice_type: str) -> float:
    """
    Calculate the lattice parameter of a crystal structure.

    Parameters
    ----------
    crystal_system : atomrdf.System
        The crystal structure as a atomrdf.System object.
    lattice_type : str
        The type of lattice ('fcc', 'bcc', 'hcp').

    Returns
    -------
    Tuple[float, float, float]
        The lattice constants of the crystal structure.

    Raises
    ------
    ValueError
        If the lattice type is not supported.
    """
    val, dist = crystal_system.calculate.radial_distribution_function(bins=500)
    peaks, _ = find_peaks(val, height=0)

    lattice_calculation = {
        "fcc": lambda d: d * sqrt(2),
        "bcc": lambda d: d * 2 / sqrt(3),
        "hcp": lambda d: d,
        "cubic diamond": lambda d: d * 4 / sqrt(3),
    }

    if lattice_type not in lattice_calculation:
        raise ValueError(f"Lattice type '{lattice_type}' not supported")

    lattice_constants = float(lattice_calculation[lattice_type](dist[peaks[0]]))
    # Round to 3 decimal places
    return round(lattice_constants, 3)

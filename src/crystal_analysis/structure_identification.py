"""Crystal structure identification and lattice parameter calculation."""

import logging
from math import sqrt

from atomrdf import System
from scipy.signal import find_peaks


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

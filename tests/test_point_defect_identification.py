"""Tests for point defect identification."""

import numpy as np
import pytest
from atomid.point_defect_analysis.wigner_seitz_method import (
    analyze_defects,
    create_index_finder,
)


@pytest.fixture
def reference_array() -> np.ndarray:
    return np.array([[0.0, 0.0], [1.0, 1.0], [2.0, 2.0], [3.0, 3.0]])


def test_analyze_defects_vacacy() -> None:
    actual_position = [
        (0.0, 0.0, 0.0),
        (2.1, 2.1, 2.1),
        (3.0, 3.0, 3.0),
        (4.0, 4.0, 4.0),
    ]
    reference_position = [
        (0.0, 0.0, 0.0),
        (1.0, 1.0, 1.0),
        (2.0, 2.0, 2.0),
        (3.0, 3.0, 3.0),
        (
            4.0,
            4.0,
            4.0,
        ),
    ]

    defects = analyze_defects(
        reference_positions=reference_position, actual_positions=actual_position
    )

    assert defects["vacancies"]["count"] == 1
    assert defects["vacancies"]["fraction"] == 0.2


def test_analyze_defects_interstitial() -> None:
    actual_position = [
        (0.0, 0.0, 0.0),
        (1.0, 1.0, 1.0),
        (1.1, 1.1, 1.1),
        (2.0, 2.0, 2.0),
        (4.0, 4.0, 4.0),
    ]
    reference_position = [
        (0.0, 0.0, 0.0),
        (1.0, 1.0, 1.0),
        (2.0, 2.0, 2.0),
        (3.0, 3.0, 3.0),
    ]

    defects = analyze_defects(
        reference_positions=reference_position, actual_positions=actual_position
    )
    assert defects["interstitials"]["count"] == 1
    assert defects["interstitials"]["fraction"] == 0.25


def test_analyze_defects_substitution() -> None:
    actual_positions = [
        (0.0, 0.0, 0.0),
        (1.0, 1.0, 1.0),
        (2.0, 2.0, 2.0),
        (3.0, 3.0, 3.0),
    ]

    reference_positions = [
        (0.0, 0.0, 0.0),
        (1.0, 1.0, 1.0),
        (2.0, 2.0, 2.0),
        (3.0, 3.0, 3.0),
    ]
    species_ref = ["H", "He", "Li", "Be"]
    species_actual = ["H", "He", "Li", "B"]

    defects = analyze_defects(
        reference_positions, actual_positions, species_ref, species_actual
    )
    assert defects["substitutions"]["count"] == 1
    assert defects["substitutions"]["fraction"] == 0.25


def test_create_index_finder_euclidean(reference_array: np.ndarray) -> None:
    finder = create_index_finder(reference_array)
    assert finder(np.array([1, 1])) == 1
    assert finder(np.array([2.1, 2.1])) == 2
    assert finder(np.array([3, 3])) == 3
    assert finder(np.array([0, 0])) == 0


def test_create_index_finder_kd_tree(reference_array: np.ndarray) -> None:
    finder = create_index_finder(reference_array, method="kd_tree")
    assert finder(np.array([1, 1])) == 1
    assert finder(np.array([2.1, 2.1])) == 2
    assert finder(np.array([3, 3])) == 3
    assert finder(np.array([0, 0])) == 0


def test_create_index_finder_annoy(reference_array: np.ndarray) -> None:
    finder = create_index_finder(reference_array, method="annoy")
    assert finder(np.array([1, 1])) == 1
    assert finder(np.array([2.1, 2.1])) == 2
    assert finder(np.array([3, 3])) == 3
    assert finder(np.array([0, 0])) == 0

"""Tests for the AnnotateCrystals class."""

from typing import List

import atomrdf as ardf
import pytest
from atomid.annotate import AnnotateCrystal
from atomid.annotate_crystals import AnnotateCrystals
from rdflib import URIRef


@pytest.fixture
def annotate_crystals() -> AnnotateCrystals:
    data_files: List[str] = [
        "tests/data/bcc/Fe/defect/interstitial/initial/Fe_interstitial.poscar",
        "tests/data/fcc/Al/no_defect/initial/Al.poscar",
    ]
    return AnnotateCrystals(data_files=data_files)


class TestAnnotateCrystals:
    """Tests for the AnnotateCrystals class."""

    def test_num_samples(self, annotate_crystals: AnnotateCrystals) -> None:
        assert annotate_crystals.num_samples == 2

    def test_kg(self, annotate_crystals: AnnotateCrystals) -> None:
        assert isinstance(annotate_crystals.kg, ardf.KnowledgeGraph)

    def test_add_data_file(self, annotate_crystals: AnnotateCrystals) -> None:
        new_file: str = "tests/data/hcp/Mg/no_defect/initial/Mg.poscar"
        annotate_crystals.add_data_file(new_file)
        assert annotate_crystals.num_samples == 3
        assert annotate_crystals.get_data_file(3) == new_file

    def test_get_crystal(self, annotate_crystals: AnnotateCrystals) -> None:
        crystal: AnnotateCrystal = annotate_crystals.get_crystal(1)
        assert isinstance(crystal, AnnotateCrystal)

    def test_get_data_file(self, annotate_crystals: AnnotateCrystals) -> None:
        data_file: str = annotate_crystals.get_data_file(1)
        assert (
            data_file
            == "tests/data/bcc/Fe/defect/interstitial/initial/Fe_interstitial.poscar"
        )

    def test_get_sample(self, annotate_crystals: AnnotateCrystals) -> None:
        sample = annotate_crystals.get_sample(1)
        assert (
            sample["data_file"]
            == "tests/data/bcc/Fe/defect/interstitial/initial/Fe_interstitial.poscar"
        )
        assert isinstance(sample["system"], AnnotateCrystal)

    def test_annotate_all_crystal_structures(
        self, annotate_crystals: AnnotateCrystals
    ) -> None:
        annotate_crystals.annotate_all_crystal_structures()
        for idx in annotate_crystals.crystals_dict:
            system: AnnotateCrystal = annotate_crystals.crystals_dict[idx]["system"]
            assert system.crystal_type is not None

    def test_annotate_crystal_structure(
        self, annotate_crystals: AnnotateCrystals
    ) -> None:
        idx: int = 1
        annotate_crystals.annotate_crystal_structure(idx)
        system: AnnotateCrystal = annotate_crystals.crystals_dict[idx]["system"]
        assert system.crystal_type is not None

    def test_identify_point_defects(self, annotate_crystals: AnnotateCrystals) -> None:
        idx = 1
        reference_file = "tests/data/bcc/Fe/no_defect/initial/Fe.poscar"
        annotate_crystals.annotate_crystal_structure(idx)
        annotate_crystals.identify_point_defects(idx, reference_file)
        system = annotate_crystals.get_crystal(idx)
        assert "interstitials" in system.defects

    def test_annotate_point_defects(self, annotate_crystals: AnnotateCrystals) -> None:
        idx = 1
        reference_file = "tests/data/bcc/Fe/no_defect/initial/Fe.poscar"
        annotate_crystals.annotate_crystal_structure(idx)
        annotate_crystals.identify_point_defects(idx, reference_file)
        annotate_crystals.annotate_point_defects(idx)

        kg = annotate_crystals.kg
        has_impurity_concentration = URIRef(
            "http://purls.helmholtz-metadaten.de/podo/hasImpurityConcentration"
        )
        has_number_of_impurity_atoms = URIRef(
            "http://purls.helmholtz-metadaten.de/podo/hasNumberOfImpurityAtoms"
        )

        impurity_concentration_exists = (
            None,
            has_impurity_concentration,
            None,
        ) in kg.graph
        number_of_impurity_atoms_exists = (
            None,
            has_number_of_impurity_atoms,
            None,
        ) in kg.graph

        assert impurity_concentration_exists
        assert number_of_impurity_atoms_exists

    def test_annotate_all_point_defects(
        self, annotate_crystals: AnnotateCrystals
    ) -> None:
        reference_files = [
            "tests/data/bcc/Fe/no_defect/initial/Fe.poscar",
            "tests/data/fcc/Al/no_defect/initial/Al.poscar",
        ]
        annotate_crystals.annotate_all_crystal_structures()
        annotate_crystals.identify_point_defects_all_samples(
            reference_files, reference_format_list=None
        )
        annotate_crystals.annotate_point_defects_all_samples()

        kg = annotate_crystals.kg
        has_impurity_concentration = URIRef(
            "http://purls.helmholtz-metadaten.de/podo/hasImpurityConcentration"
        )
        has_number_of_impurity_atoms = URIRef(
            "http://purls.helmholtz-metadaten.de/podo/hasNumberOfImpurityAtoms"
        )

        impurity_concentration_exists = (
            None,
            has_impurity_concentration,
            None,
        ) in kg.graph
        number_of_impurity_atoms_exists = (
            None,
            has_number_of_impurity_atoms,
            None,
        ) in kg.graph

        assert impurity_concentration_exists
        assert number_of_impurity_atoms_exists

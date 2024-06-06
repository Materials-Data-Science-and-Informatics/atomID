"""Tests for the AnnotateCrystal class in atomid/annotate.py."""

import os

import pytest
from atomid.annotate import AnnotateCrystal
from utils.compare_rdf import compare_graphs


class TestAnnotateCrystal:
    """Tests for the AnnotateCrystal class."""

    testing_dict = {
        "fcc": "Al",
        "bcc": "Fe",
        "hcp": "Mg",
        "diamond": "Si",
    }

    # Generate combinations of sample and reference crystal files
    test_data_combinations = []

    defect_types = ["interstitial", "substitution", "vacancy"]

    for structure, element in testing_dict.items():
        for defect in defect_types:
            test_data_combinations.append(
                (
                    f"tests/data/{structure}/{element}/defect/{defect}/initial/{element}_{defect}.poscar",
                    f"tests/data/{structure}/{element}/no_defect/initial/{element}.poscar",
                )
            )

    @pytest.mark.parametrize(
        "sample_crystal_file, reference_crystal_file", test_data_combinations
    )
    def test_read_crystal_structure_file(
        self, sample_crystal_file: str, reference_crystal_file: str
    ) -> None:
        annotate_crystal = AnnotateCrystal()
        annotate_crystal.read_crystal_structure_file(
            sample_crystal_file, format="vasp"
        )  # Adjust format if needed

        assert annotate_crystal.ase_crystal is not None
        assert annotate_crystal.system is not None
        assert annotate_crystal.kg is not None

    @pytest.mark.parametrize(
        "sample_crystal_file, reference_crystal_file", test_data_combinations
    )
    def test_annotate_crystal_structure(
        self, sample_crystal_file: str, reference_crystal_file: str
    ) -> None:
        annotate_crystal = AnnotateCrystal()
        annotate_crystal.read_crystal_structure_file(
            sample_crystal_file, format="vasp"
        )  # Adjust format if needed
        annotate_crystal.annotate_crystal_structure()

        assert annotate_crystal.system is not None

    @pytest.mark.parametrize(
        "sample_crystal_file, reference_crystal_file", test_data_combinations
    )
    def test_identify_defects(
        self, sample_crystal_file: str, reference_crystal_file: str
    ) -> None:
        annotate_crystal = AnnotateCrystal()
        annotate_crystal.read_crystal_structure_file(
            sample_crystal_file, format="vasp"
        )  # Adjust format if needed
        defects = annotate_crystal.identify_point_defects(
            reference_crystal_file, ref_format="vasp"
        )  # Adjust format if needed

        assert isinstance(defects, dict)
        assert "vacancies" in defects
        assert "interstitials" in defects
        assert "substitutions" in defects

    @pytest.mark.parametrize(
        "sample_crystal_file, reference_crystal_file", test_data_combinations
    )
    def test_write_defects(
        self, tmp_path: str, sample_crystal_file: str, reference_crystal_file: str
    ) -> None:
        """Test writing the defects to a file."""
        annotate_crystal = AnnotateCrystal()
        annotate_crystal.read_crystal_structure_file(sample_crystal_file, format="vasp")
        annotate_crystal.identify_point_defects(reference_crystal_file, ref_format="vasp")

        annotate_crystal.write_to_file(f"{tmp_path}/annotated_output.ttl", "ttl")

        assert os.path.exists(f"{tmp_path}/annotated_output.ttl")

    @pytest.mark.parametrize(
        "sample_crystal_file, reference_crystal_file", test_data_combinations
    )
    def test_output_annotation(
        self, tmp_path: str, sample_crystal_file: str, reference_crystal_file: str
    ) -> None:
        """Test the output of the annotation."""
        annotate_crystal = AnnotateCrystal()
        annotate_crystal.read_crystal_structure_file(sample_crystal_file, format="vasp")
        annotate_crystal.annotate_point_defects(reference_crystal_file, ref_format="vasp")
        annotate_crystal.write_to_file(f"{tmp_path}/annotated_output.ttl", "ttl")

        assert os.path.exists(f"{tmp_path}/annotated_output.ttl")

        result, differences = compare_graphs(
            f"{tmp_path}/annotated_output.ttl",
            sample_crystal_file.replace("poscar", "ttl"),
        )

        if result is False:
            print(differences)

        assert result


if __name__ == "__main__":
    pytest.main()

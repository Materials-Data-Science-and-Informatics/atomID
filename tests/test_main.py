"""Tests for the AnnotateCrystal class in atomid/annotate.py."""

import os

import pytest
from atomid.annotate import AnnotateCrystal
from utils.compare_rdf import compare_graphs
import tempfile



class TestAnnotateCrystal:
    """Tests for the AnnotateCrystal class."""

    # Define the combinations of sample and reference crystal files
    test_data_combinations = [
        (
            "tests/data/fcc/Al/defect/interstitial/initial/Al_interstitial.poscar",
            "tests/data/fcc/Al/no_defect/initial/Al.poscar",
        ),
    ]
    testing_dict = {
        "fcc": "Al",
        "bcc": "Fe",
        "hcp": "Mg",
        "diamond": "Si",
    }

    @pytest.mark.parametrize(
        "sample_crystal_file, reference_crystal_file", test_data_combinations
    )
    def test_read_crystal_structure_file(
        self, sample_crystal_file, reference_crystal_file
    ):
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
        self, sample_crystal_file, reference_crystal_file
    ):
        annotate_crystal = AnnotateCrystal()
        annotate_crystal.read_crystal_structure_file(
            sample_crystal_file, format="vasp"
        )  # Adjust format if needed
        annotate_crystal.annotate_crystal_structure()

        assert annotate_crystal.system is not None

    @pytest.mark.parametrize(
        "sample_crystal_file, reference_crystal_file", test_data_combinations
    )
    def test_identify_defects(self, sample_crystal_file, reference_crystal_file):
        annotate_crystal = AnnotateCrystal()
        annotate_crystal.read_crystal_structure_file(
            sample_crystal_file, format="vasp"
        )  # Adjust format if needed
        defects = annotate_crystal.identify_defects(
            reference_crystal_file, ref_format="vasp"
        )  # Adjust format if needed

        assert isinstance(defects, dict)
        assert "Vacancies" in defects
        assert "Interstitials" in defects
        assert "Substitutions" in defects

    @pytest.mark.parametrize(
        "sample_crystal_file, reference_crystal_file", test_data_combinations
    )
    def test_write_defects(self, sample_crystal_file, reference_crystal_file):
        """Test writing the defects to a file."""
        annotate_crystal = AnnotateCrystal()
        annotate_crystal.read_crystal_structure_file(sample_crystal_file, format="vasp")
        annotate_crystal.identify_defects(reference_crystal_file, ref_format="vasp")

        tmpdir = tempfile.mkdtemp()
        annotate_crystal.write_to_file(f"{tmpdir}/annotated_output.ttl", "ttl")

        assert os.path.exists(f"{tmpdir}/annotated_output.ttl")

        # Clean up
        os.remove(f"{tmpdir}/annotated_output.ttl")

    @pytest.mark.parametrize(
        "sample_crystal_file, reference_crystal_file", test_data_combinations
    )
    def test_output_annotation(self, sample_crystal_file, reference_crystal_file):
        """Test the output of the annotation."""
        annotate_crystal = AnnotateCrystal()
        annotate_crystal.read_crystal_structure_file(sample_crystal_file, format="vasp")
        annotate_crystal.identify_defects(
            reference_crystal_file, ref_format="vasp"
        )

        tmpdir = tempfile.mkdtemp()
        annotate_crystal.write_to_file(f"{tmpdir}/annotated_output.ttl", "ttl")

        assert os.path.exists(f"{tmpdir}/annotated_output.ttl")

        result, differences = compare_graphs(
            f"{tmpdir}/annotated_output.ttl",
            "tests/data/fcc/Al/defect/interstitial/initial/Al_interstitial.ttl",
        )

        if result is False:
            print(differences)

        assert result
        # Clean up
        os.remove(f"{tmpdir}/annotated_output.ttl")


if __name__ == "__main__":
    pytest.main()

"""Tests for the AnnotateCrystal class in atomid/annotate.py."""

import os

import pytest
from atomid.annotate import AnnotateCrystal
from utils.compare_rdf import compare_graphs


class TestAnnotateSystemTest:
    """Tests for the AnnotateCrystal class for Point defects."""

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
        annotate_crystal.identify_crystal_structure()
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
        annotate_crystal.identify_point_defects(
            reference_crystal_file, ref_format="vasp"
        )
        annotate_crystal.annotate_point_defects()
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

        annotate_crystal.identify_crystal_structure()
        annotate_crystal.annotate_crystal_structure()

        annotate_crystal.identify_point_defects(
            reference_crystal_file, ref_format="vasp"
        )
        annotate_crystal.annotate_point_defects()
        annotate_crystal.write_to_file(f"{tmp_path}/annotated_output.ttl", "ttl")

        assert os.path.exists(f"{tmp_path}/annotated_output.ttl")

        result, differences = compare_graphs(
            f"{tmp_path}/annotated_output.ttl",
            sample_crystal_file.replace("poscar", "ttl"),
        )

        if result is False:
            print(differences)

        assert result


class TestAnnotatePointDefects:
    """Tests for the AnnotateCrystal class for Point defects."""

    def test_set_lattice_constant(self) -> None:
        annotate_crystal = AnnotateCrystal(
            "tests/data/fcc/Al/no_defect/initial/Al.poscar", "vasp"
        )
        annotate_crystal.set_lattice_constant(3.5)
        assert annotate_crystal.lattice_constant == 3.5

    def test_set_crystal_structure(self) -> None:
        annotate_crystal = AnnotateCrystal(
            "tests/data/fcc/Al/no_defect/initial/Al.poscar", "vasp"
        )
        annotate_crystal.set_crystal_structure("fcc")
        assert annotate_crystal.crystal_type == "fcc"

    def test_annotate_returns_error_no_lattice_constant(self) -> None:
        annotate_crystal = AnnotateCrystal()
        annotate_crystal.read_crystal_structure_file(
            "tests/data/fcc/Al/no_defect/initial/Al.poscar", format="vasp"
        )
        annotate_crystal.set_crystal_structure("fcc")

        with pytest.raises(ValueError):
            annotate_crystal.annotate_crystal_structure()

    def test_annotate_returns_error_no_crystal_type(self) -> None:
        annotate_crystal = AnnotateCrystal()
        annotate_crystal.read_crystal_structure_file(
            "tests/data/fcc/Al/no_defect/initial/Al.poscar", format="vasp"
        )
        annotate_crystal.set_lattice_constant(3.5)
        with pytest.raises(ValueError):
            annotate_crystal.annotate_crystal_structure()

    def test_add_vacancy(self) -> None:
        annotate_crystal = AnnotateCrystal()
        annotate_crystal.read_crystal_structure_file(
            "tests/data/fcc/Al/defect/vacancy/initial/Al_vacancy.poscar", format="vasp"
        )
        annotate_crystal.identify_crystal_structure()
        annotate_crystal.annotate_crystal_structure()
        annotate_crystal.add_vacancy_information(concentration=0.1, number=40)
        assert annotate_crystal.defects["vacancies"]["count"] == 40
        assert annotate_crystal.defects["vacancies"]["concentration"] == 0.1

    def test_add_interstitial(self) -> None:
        annotate_crystal = AnnotateCrystal()
        annotate_crystal.read_crystal_structure_file(
            "tests/data/fcc/Al/defect/interstitial/initial/Al_interstitial.poscar",
            format="vasp",
        )
        annotate_crystal.identify_crystal_structure()
        annotate_crystal.annotate_crystal_structure()
        annotate_crystal.add_interstitial_information(concentration=0.1, number=40)
        assert annotate_crystal.defects["interstitials"]["count"] == 40
        assert annotate_crystal.defects["interstitials"]["concentration"] == 0.1

    def test_add_substitution(self) -> None:
        annotate_crystal = AnnotateCrystal()
        annotate_crystal.read_crystal_structure_file(
            "tests/data/fcc/Al/defect/substitution/initial/Al_substitution.poscar",
            format="vasp",
        )
        annotate_crystal.identify_crystal_structure()
        annotate_crystal.annotate_crystal_structure()
        annotate_crystal.add_substitution_information(concentration=0.1, number=40)
        assert annotate_crystal.defects["substitutions"]["count"] == 40
        assert annotate_crystal.defects["substitutions"]["concentration"] == 0.1


class TestAnnotateGrains:
    """Tests for the AnnotateCrystal class for Grain boundaries."""

    def test_identify_grain_boundary(self) -> None:
        grain_boundary_path = "tests/data/grain_boundary/fe_grain_boundary.poscar"

        annotate_crystal = AnnotateCrystal()
        annotate_crystal.read_crystal_structure_file(grain_boundary_path, format="vasp")
        annotate_crystal.identify_crystal_structure()
        annotate_crystal.annotate_crystal_structure()
        grains, angles = annotate_crystal.identify_grains()

        assert grains is not None
        assert angles is not None


class TestAnnotateDislocations:
    """Tests for the AnnotateCrystal class for Line defects."""

    def test_identify_dislocations(self) -> None:
        dislocation_path = "tests/data/dislocation/Al_edge.cfg"

        annotate_crystal = AnnotateCrystal()
        annotate_crystal.read_crystal_structure_file(dislocation_path, format="cfg")
        annotate_crystal.identify_crystal_structure()
        annotate_crystal.annotate_crystal_structure()
        burgers_vectors, lengths = annotate_crystal.identify_line_defects()

        assert burgers_vectors is not None
        assert lengths is not None

    def test_identify_dislocations_none(self) -> None:
        dislocation_path = "tests/data/fcc/Al/no_defect/initial/Al.poscar"

        annotate_crystal = AnnotateCrystal()
        annotate_crystal.read_crystal_structure_file(dislocation_path, format="vasp")
        annotate_crystal.identify_crystal_structure()
        annotate_crystal.annotate_crystal_structure()
        burgers_vectors, lengths = annotate_crystal.identify_line_defects()

        assert burgers_vectors is None
        assert lengths is None


if __name__ == "__main__":
    pytest.main()

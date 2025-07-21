# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scienting Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

import pytest

from .data_validator import DataValidator
from .synthetic_data import SyntheticAtomDataGenerator


@pytest.fixture
def temp_project_dir(tmp_path):
    """Provides a temporary directory for project stores."""
    dir_path = tmp_path / "atomea_synthetic_test_project"
    dir_path.mkdir()
    yield str(dir_path)
    # Clean up is handled by tmp_path fixture itself


def test_generate_and_validate_valid_data(temp_project_dir):
    """
    Test generating standard, valid synthetic data and validating it.
    """
    generator = SyntheticAtomDataGenerator(
        n_frames=(5, 10),
        n_atoms=(20, 50),
        random_seed=42,  # Ensure reproducibility for this test
    )

    # Generate and write data to a new project
    project, reference_data = generator.create_synthetic_project(
        temp_project_dir,
        ens_id="test_ensemble",
        run_id="run_valid",
        clean_up_paths=True,
    )

    # Initialize validator and run all checks
    validator = DataValidator(project, reference_data, expect_fuzzed_data=False)
    validator.validate_all_data()


def test_fuzz_invalid_coordinates_shape(temp_project_dir):
    """
    Test fuzzing with invalid coordinate shapes.
    Expects an error during data writing or reading due to shape mismatch.
    """
    generator = SyntheticAtomDataGenerator(
        n_frames=1,
        n_atoms=10,
        inject_invalid_shape_prob=1.0,  # Always inject invalid shape
        random_seed=43,
    )

    # We expect an error when writing or reading due to the malformed shape.
    # The actual point of failure depends on the underlying store/Data object.
    # We wrap the call in pytest.raises to capture the expected error.
    with pytest.raises(Exception) as excinfo:
        project, reference_data = generator.create_synthetic_project(
            temp_project_dir,
            ens_id="test_fuzz_coords_shape",
            run_id="run_fuzz_coords_shape",
            clean_up_paths=True,
        )
        # Attempt to read to trigger any lazy shape validation in Zarr/Numpy
        # Note: If your Zarr/Data object doesn't immediately validate shape on write/read,
        # the error might only appear when DataValidator attempts to check `coords.shape`.
        _ = project.get_ensemble("test_fuzz_coords_shape").coordinates.read(
            run_id="run_fuzz_coords_shape"
        )

        # If the above didn't raise, try to validate, which should definitely fail on shape.
        validator = DataValidator(project, reference_data, expect_fuzzed_data=True)
        validator.validate_coordinates()  # This should fail on shape

    # You can be more specific with the exception type if your 'atomea' library
    # raises specific errors for shape mismatches, e.g., AssertionError, ValueError, etc.
    print(
        f"\nCaught expected error during fuzzing invalid coordinate shape: {excinfo.type.__name__}: {excinfo.value}"
    )


def test_fuzz_nan_coordinates(temp_project_dir):
    """
    Test fuzzing with NaN values in coordinates.
    The validator should detect NaNs unless explicitly told to ignore.
    """
    generator = SyntheticAtomDataGenerator(
        n_frames=5,
        n_atoms=10,
        inject_nan_coords_prob=1.0,  # Always inject NaNs
        random_seed=44,
    )

    project, reference_data = generator.create_synthetic_project(
        temp_project_dir,
        ens_id="test_fuzz_nan_coords",
        run_id="run_fuzz_nan_coords",
        clean_up_paths=True,
    )

    # When expect_fuzzed_data is False, validator should detect NaNs and raise.
    with pytest.raises(AssertionError, match="NaN values found in Coordinates."):
        validator_strict = DataValidator(
            project, reference_data, expect_fuzzed_data=False
        )
        validator_strict.validate_coordinates()

    print("\nCaught expected AssertionError for NaN coordinates in strict validation.")

    # When expect_fuzzed_data is True, validator should not raise for NaNs (as configured).
    # The _validate_array_data helper in DataValidator has `check_nan_inf=not self.expect_fuzzed_data`
    # so this test case primarily validates that the *fuzzing parameter* influences the validator.
    validator_fuzz = DataValidator(project, reference_data, expect_fuzzed_data=True)
    # This call should pass because check_nan_inf will be False
    validator_fuzz.validate_coordinates()
    print(
        "Coordinates with NaNs passed validation when expect_fuzzed_data is True (NaN check skipped)."
    )


def test_fuzz_out_of_bounds_atom_id(temp_project_dir):
    """
    Test fuzzing with out-of-bounds atomic numbers (uint8 overflow).
    Expects an error or incorrect data when reading, depending on NumPy's behavior.
    """
    generator = SyntheticAtomDataGenerator(
        n_atoms=5,
        inject_out_of_bounds_id_prob=1.0,  # Always inject out-of-bounds ID
        random_seed=45,
    )

    with pytest.raises(
        OverflowError, match=r"Python integer \d+ out of bounds for uint8"
    ):
        project, reference_data = generator.create_synthetic_project(
            temp_project_dir,
            ens_id="test_fuzz_id_bounds",
            run_id="run_fuzz_id_bounds",
            clean_up_paths=True,
        )

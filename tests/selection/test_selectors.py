# tests/selection/test_selectors.py
import logging
from pathlib import Path

import numpy as np
import pytest

from atomea.containers import Ensemble, Project
from atomea.selection.expressions import (
    AndExpression,
    AtomTypeIs,
    DistanceWithin,
    MolIdIs,
    NotExpression,
    OrExpression,
)
from atomea.selection.selector import EnsembleSelector
from atomea.stores import DiskFormat
from atomea.stores.arrays import ZarrArrayStore
from atomea.stores.tables import PolarsTableStore

# Configure logging for capturing warnings/debug messages
logger = logging.getLogger("atomea")


def assert_masks_equal(actual_mask_iterator, expected_masks, msg=""):
    """Helper to compare iterators of masks."""
    for i, (actual, expected) in enumerate(zip(actual_mask_iterator, expected_masks)):
        np.testing.assert_array_equal(
            actual, expected, f"Mask mismatch at microstate {i}. {msg}"
        )
    # Ensure both iterators are exhausted
    with pytest.raises(StopIteration):
        next(actual_mask_iterator)


class TestSelectionExpressions:
    """Tests for individual SelectionExpression classes."""

    def test_atom_types_single_type(
        self, test_project_with_synthetic_data: tuple[Project, Ensemble, str, int, int]
    ):
        project, ensemble, run_id, n_atoms, n_frames = test_project_with_synthetic_data

        # Select Oxygen atoms ("OW" type)
        # Atoms: O0, H1, H2, O3, H4, H5, C6, H7, H8, H9, H10
        # Types: OW, HW, HW, OW, HW, HW, C, H, H, H, H
        selector = AtomTypeIs(atom_types=["OW"])
        expected_mask = np.array(
            [True, False, False, True, False, False, False, False, False, False, False],
            dtype=np.bool_,
        )
        # AtomTypeIs is ENSEMBLE cadence, so mask is same for all frames
        expected_masks = [expected_mask] * n_frames

        actual_masks_iterator = selector.evaluate(
            ensemble=ensemble, run_id=run_id, micro_id=slice(None)
        )
        assert_masks_equal(actual_masks_iterator, expected_masks)

        # Test single microstate
        actual_mask_iterator_single = selector.evaluate(
            ensemble=ensemble, run_id=run_id, micro_id=0
        )
        assert_masks_equal(actual_mask_iterator_single, [expected_mask])

    def test_atom_types_multiple_types(
        self, test_project_with_synthetic_data: tuple[Project, Ensemble, str, int, int]
    ):
        project, ensemble, run_id, n_atoms, n_frames = test_project_with_synthetic_data

        # Select all Hydrogen atoms ("HW" and "H" types)
        # Atoms: O0, H1, H2, O3, H4, H5, C6, H7, H8, H9, H10
        # Types: OW, HW, HW, OW, HW, HW, C, H, H, H, H
        selector = AtomTypeIs(atom_types=["HW", "H"])
        expected_mask = np.array(
            [False, True, True, False, True, True, False, True, True, True, True],
            dtype=np.bool_,
        )
        expected_masks = [expected_mask] * n_frames

        actual_masks_iterator = selector.evaluate(
            ensemble=ensemble, run_id=run_id, micro_id=slice(None)
        )
        assert_masks_equal(actual_masks_iterator, expected_masks)

    def test_atom_types_no_match(
        self, test_project_with_synthetic_data: tuple[Project, Ensemble, str, int, int]
    ):
        project, ensemble, run_id, n_atoms, n_frames = test_project_with_synthetic_data

        # Select a non-existent atom type
        selector = AtomTypeIs(atom_types=["XYZ"])
        expected_mask = np.zeros(n_atoms, dtype=np.bool_)
        expected_masks = [expected_mask] * n_frames

        actual_masks_iterator = selector.evaluate(
            ensemble=ensemble, run_id=run_id, micro_id=slice(None)
        )
        assert_masks_equal(actual_masks_iterator, expected_masks)

    def test_atom_types_empty_list(
        self, test_project_with_synthetic_data: tuple[Project, Ensemble, str, int, int]
    ):
        project, ensemble, run_id, n_atoms, n_frames = test_project_with_synthetic_data

        # Select with an empty list of atom types
        selector = AtomTypeIs(atom_types=[])
        expected_mask = np.zeros(n_atoms, dtype=np.bool_)
        expected_masks = [expected_mask] * n_frames

        actual_masks_iterator = selector.evaluate(
            ensemble=ensemble, run_id=run_id, micro_id=slice(None)
        )
        assert_masks_equal(actual_masks_iterator, expected_masks)

    def test_molecule_ids_single_id(
        self, test_project_with_synthetic_data: tuple[Project, Ensemble, str, int, int]
    ):
        project, ensemble, run_id, n_atoms, n_frames = test_project_with_synthetic_data

        # Select molecule ID 0 (first water)
        # Mol IDs: 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 2
        selector = MolIdIs(mol_ids=[0])
        expected_mask = np.array(
            [True, True, True, False, False, False, False, False, False, False, False],
            dtype=np.bool_,
        )
        expected_masks = [expected_mask] * n_frames

        actual_masks_iterator = selector.evaluate(
            ensemble=ensemble, run_id=run_id, micro_id=slice(None)
        )
        assert_masks_equal(actual_masks_iterator, expected_masks)

    def test_molecule_ids_multiple_ids(
        self, test_project_with_synthetic_data: tuple[Project, Ensemble, str, int, int]
    ):
        project, ensemble, run_id, n_atoms, n_frames = test_project_with_synthetic_data

        # Select molecule IDs 0 and 2 (first water and methane)
        selector = MolIdIs(mol_ids=[0, 2])
        expected_mask = np.array(
            [True, True, True, False, False, False, True, True, True, True, True],
            dtype=np.bool_,
        )
        expected_masks = [expected_mask] * n_frames

        actual_masks_iterator = selector.evaluate(
            ensemble=ensemble, run_id=run_id, micro_id=slice(None)
        )
        assert_masks_equal(actual_masks_iterator, expected_masks)

    def test_molecule_ids_no_match(
        self, test_project_with_synthetic_data: tuple[Project, Ensemble, str, int, int]
    ):
        project, ensemble, run_id, n_atoms, n_frames = test_project_with_synthetic_data

        # Select a non-existent molecule ID
        selector = MolIdIs(mol_ids=[99])
        expected_mask = np.zeros(n_atoms, dtype=np.bool_)
        expected_masks = [expected_mask] * n_frames

        actual_masks_iterator = selector.evaluate(
            ensemble=ensemble, run_id=run_id, micro_id=slice(None)
        )
        assert_masks_equal(actual_masks_iterator, expected_masks)

    def test_distance_within_from_slice(
        self, test_project_with_synthetic_data: tuple[Project, Ensemble, str, int, int]
    ):
        project, ensemble, run_id, n_atoms, n_frames = test_project_with_synthetic_data

        # Select atoms within 1.0 A of atom 0 (O of first water)
        # Atom 0: O at (0.0, 0.0, 0.0)
        # Atom 1: H at (0.7, 0.7, 0.0) -> dist = sqrt(0.7^2 + 0.7^2) = sqrt(0.98) approx 0.99
        # Atom 2: H at (0.7, -0.7, 0.0) -> dist = sqrt(0.7^2 + (-0.7)^2) = sqrt(0.98) approx 0.99
        # Both H1 and H2 are within 1.0 A of O0. O0 is also selected.
        # Other atoms are far away.
        selector = DistanceWithin(
            from_atoms=slice(0, 1), dist=1.0
        )  # Atom 0 is the reference
        expected_mask_ms0 = np.array(
            [True, True, True, False, False, False, False, False, False, False, False],
            dtype=np.bool_,
        )
        # For MS1, atom 0 is at (0.1, 0.0, 0.0)
        # Atom 1 is at (0.8, 0.7, 0.0). dist = sqrt((0.8-0.1)^2 + 0.7^2) = sqrt(0.7^2 + 0.7^2) = 0.99
        # Atom 2 is at (0.8, -0.7, 0.0). dist = 0.99
        expected_mask_ms1 = np.array(
            [True, True, True, False, False, False, False, False, False, False, False],
            dtype=np.bool_,
        )
        expected_masks = [expected_mask_ms0, expected_mask_ms1]

        actual_masks_iterator = selector.evaluate(
            ensemble=ensemble, run_id=run_id, micro_id=slice(None)
        )
        assert_masks_equal(actual_masks_iterator, expected_masks)

        # Test single microstate
        actual_mask_iterator_single = selector.evaluate(
            ensemble=ensemble, run_id=run_id, micro_id=0
        )
        assert_masks_equal(actual_mask_iterator_single, [expected_mask_ms0])

    def test_distance_within_from_selection_expression(
        self, test_project_with_synthetic_data: tuple[Project, Ensemble, str, int, int]
    ):
        project, ensemble, run_id, n_atoms, n_frames = test_project_with_synthetic_data

        # Select atoms within 1.0 A of any "OW" (Oxygen from Water) atoms
        # "OW" atoms are indices 0 and 3.
        # MS0:
        # Ref 0: O0 (0,0,0) -> H1, H2 are within 1.0
        # Ref 3: O3 (3,0,0) -> H4, H5 are within 1.0
        # Expected selection: [O0, H1, H2, O3, H4, H5]
        ref_selection = AtomTypeIs(atom_types=["OW"])
        selector = DistanceWithin(from_atoms=ref_selection, dist=1.0)
        expected_mask = np.array(
            [True, True, True, True, True, True, False, False, False, False, False],
            dtype=np.bool_,
        )
        expected_masks = [
            expected_mask,
            expected_mask,
        ]  # Shifted positions maintain same relative distances

        actual_masks_iterator = selector.evaluate(
            ensemble=ensemble, run_id=run_id, micro_id=slice(None)
        )
        assert_masks_equal(actual_masks_iterator, expected_masks)

    def test_distance_within_dist_zero(
        self, test_project_with_synthetic_data: tuple[Project, Ensemble, str, int, int]
    ):
        project, ensemble, run_id, n_atoms, n_frames = test_project_with_synthetic_data

        # Only atoms exactly at the reference point should be selected
        # Selecting around atom 0 (O0)
        selector = DistanceWithin(from_atoms=slice(0, 1), dist=0.0)
        expected_mask = np.array(
            [
                True,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
            ],
            dtype=np.bool_,
        )
        expected_masks = [expected_mask] * n_frames

        actual_masks_iterator = selector.evaluate(
            ensemble=ensemble, run_id=run_id, micro_id=slice(None)
        )
        assert_masks_equal(actual_masks_iterator, expected_masks)

    def test_distance_within_dist_negative(
        self, test_project_with_synthetic_data: tuple[Project, Ensemble, str, int, int]
    ):
        project, ensemble, run_id, n_atoms, n_frames = test_project_with_synthetic_data

        # Negative distance should always result in an empty selection
        selector = DistanceWithin(from_atoms=slice(0, 1), dist=-1.0)
        expected_mask = np.zeros(n_atoms, dtype=np.bool_)
        expected_masks = [expected_mask] * n_frames

        actual_masks_iterator = selector.evaluate(
            ensemble=ensemble, run_id=run_id, micro_id=slice(None)
        )
        assert_masks_equal(actual_masks_iterator, expected_masks)

    def test_distance_within_no_reference_atoms(
        self,
        test_project_with_synthetic_data: tuple[Project, Ensemble, str, int, int],
        caplog,
    ):
        project, ensemble, run_id, n_atoms, n_frames = test_project_with_synthetic_data
        caplog.set_level(logging.WARNING)

        # Reference selection yields no atoms
        no_atoms_selector = AtomTypeIs(atom_types=["NONEXISTENT"])
        selector = DistanceWithin(from_atoms=no_atoms_selector, dist=1.0)

        expected_mask = np.zeros(n_atoms, dtype=np.bool_)
        expected_masks = [expected_mask] * n_frames

        actual_masks_iterator = selector.evaluate(
            ensemble=ensemble, run_id=run_id, micro_id=slice(None)
        )
        assert_masks_equal(actual_masks_iterator, expected_masks)

        # Check if a warning was logged (depends on internal implementation)
        # As per the `DistanceWithin` snippet, it handles `from_coords.shape[0] == 0`
        # but doesn't necessarily log a warning itself. The warning might come from
        # a sub-expression if its data is missing. For this specific case,
        # no warning is expected from DistanceWithin if the *reference selection* is empty.
        # If the *coordinates data itself* is missing, then a warning might be logged.
        # Let's adjust based on the current implementation for `MolIdIs` and `AtomTypeIs`
        # which log warnings if their `_cached_all_mol_ids` or `_cached_all_atom_types` are None.
        # DistanceWithin doesn't have such a cache.

    def test_and_expression(
        self, test_project_with_synthetic_data: tuple[Project, Ensemble, str, int, int]
    ):
        project, ensemble, run_id, n_atoms, n_frames = test_project_with_synthetic_data

        # Select atoms that are Hydrogen ("HW" or "H") AND are in Molecule 0
        expr1 = AtomTypeIs(
            atom_types=["HW", "H"]
        )  # All H atoms: [F,T,T,F,T,T,F,T,T,T,T]
        expr2 = MolIdIs(mol_ids=[0])  # Mol 0 atoms: [T,T,T,F,F,F,F,F,F,F,F]

        # Expected: intersection of H atoms and Mol 0 (H1, H2)
        expected_mask = np.array(
            [False, True, True, False, False, False, False, False, False, False, False],
            dtype=np.bool_,
        )
        expected_masks = [expected_mask] * n_frames

        selector = AndExpression(expr1, expr2)
        actual_masks_iterator = selector.evaluate(
            ensemble=ensemble, run_id=run_id, micro_id=slice(None)
        )
        assert_masks_equal(actual_masks_iterator, expected_masks)

    def test_or_expression(
        self, test_project_with_synthetic_data: tuple[Project, Ensemble, str, int, int]
    ):
        project, ensemble, run_id, n_atoms, n_frames = test_project_with_synthetic_data

        # Select atoms that are Oxygen ("OW") OR are in Molecule 2 (Methane)
        expr1 = AtomTypeIs(atom_types=["OW"])  # O0, O3: [T,F,F,T,F,F,F,F,F,F,F]
        expr2 = MolIdIs(mol_ids=[2])  # Mol 2: [F,F,F,F,F,F,T,T,T,T,T]

        # Expected: union of O atoms and Mol 2 (O0, O3, C6, H7, H8, H9, H10)
        expected_mask = np.array(
            [True, False, False, True, False, False, True, True, True, True, True],
            dtype=np.bool_,
        )
        expected_masks = [expected_mask] * n_frames

        selector = OrExpression(expr1, expr2)
        actual_masks_iterator = selector.evaluate(
            ensemble=ensemble, run_id=run_id, micro_id=slice(None)
        )
        assert_masks_equal(actual_masks_iterator, expected_masks)

    def test_not_expression(
        self, test_project_with_synthetic_data: tuple[Project, Ensemble, str, int, int]
    ):
        project, ensemble, run_id, n_atoms, n_frames = test_project_with_synthetic_data

        # Select atoms that are NOT "HW" (Hydrogen in Water)
        expr = AtomTypeIs(atom_types=["HW"])  # HW atoms: [F,T,T,F,T,T,F,F,F,F,F]
        # Expected: all non-HW atoms (O0, O3, C6, H7, H8, H9, H10)
        expected_mask = np.array(
            [True, False, False, True, False, False, True, True, True, True, True],
            dtype=np.bool_,
        )
        expected_masks = [expected_mask] * n_frames

        selector = NotExpression(expr)
        actual_masks_iterator = selector.evaluate(
            ensemble=ensemble, run_id=run_id, micro_id=slice(None)
        )
        assert_masks_equal(actual_masks_iterator, expected_masks)

    def test_complex_logical_expression(
        self, test_project_with_synthetic_data: tuple[Project, Ensemble, str, int, int]
    ):
        project, ensemble, run_id, n_atoms, n_frames = test_project_with_synthetic_data

        # Select atoms that are (Hydrogen AND in Mol 0) OR (Carbon AND within 1.0 A of atom 6 (C6))
        # H in Mol 0: H1, H2
        # C in Mol 2 (C6) within 1.0 A of C6: only C6 itself
        # Expected: H1, H2, C6
        h_in_mol0 = AtomTypeIs(atom_types=["H", "HW"]) & MolIdIs(mol_ids=[0])
        carbon_near_c6 = AtomTypeIs(atom_types=["C"]) & DistanceWithin(
            from_atoms=slice(6, 7), dist=1.0
        )

        selector = h_in_mol0 | carbon_near_c6

        expected_mask = np.array(
            [False, True, True, False, False, False, True, False, False, False, False],
            dtype=np.bool_,
        )
        expected_masks = [expected_mask] * n_frames

        actual_masks_iterator = selector.evaluate(
            ensemble=ensemble, run_id=run_id, micro_id=slice(None)
        )
        assert_masks_equal(actual_masks_iterator, expected_masks)


class TestEnsembleSelectorFluentAPI:
    """Tests for the EnsembleSelector's fluent API."""

    def test_select_by_atom_type(
        self, test_project_with_synthetic_data: tuple[Project, Ensemble, str, int, int]
    ):
        project, ensemble, run_id, n_atoms, n_frames = test_project_with_synthetic_data

        selector = EnsembleSelector(ensemble, run_id).atom_types(["OW"])
        expected_mask = np.array(
            [True, False, False, True, False, False, False, False, False, False, False],
            dtype=np.bool_,
        )
        expected_masks = [expected_mask] * n_frames

        actual_masks_iterator = selector.get_mask(micro_id=slice(None))
        assert_masks_equal(actual_masks_iterator, expected_masks)

    def test_select_by_mol_id(
        self, test_project_with_synthetic_data: tuple[Project, Ensemble, str, int, int]
    ):
        project, ensemble, run_id, n_atoms, n_frames = test_project_with_synthetic_data

        selector = EnsembleSelector(ensemble, run_id).molecule_ids(
            [1]
        )  # Select second water
        expected_mask = np.array(
            [False, False, False, True, True, True, False, False, False, False, False],
            dtype=np.bool_,
        )
        expected_masks = [expected_mask] * n_frames

        actual_masks_iterator = selector.get_mask(micro_id=slice(None))
        assert_masks_equal(actual_masks_iterator, expected_masks)

    def test_select_distance_within(
        self, test_project_with_synthetic_data: tuple[Project, Ensemble, str, int, int]
    ):
        project, ensemble, run_id, n_atoms, n_frames = test_project_with_synthetic_data

        # Select atoms within 1.0 A of atom 3 (O of second water)
        selector = EnsembleSelector(ensemble, run_id).distance_within(
            from_atoms=3, dist=1.0
        )
        expected_mask = np.array(
            [False, False, False, True, True, True, False, False, False, False, False],
            dtype=np.bool_,
        )
        expected_masks = [
            expected_mask
        ] * n_frames  # Same relative distances for both microstates

        actual_masks_iterator = selector.get_mask(micro_id=slice(None))
        assert_masks_equal(actual_masks_iterator, expected_masks)

    def test_fluent_and_operation(
        self, test_project_with_synthetic_data: tuple[Project, Ensemble, str, int, int]
    ):
        project, ensemble, run_id, n_atoms, n_frames = test_project_with_synthetic_data

        # Select Carbon atoms ("C" type) AND are in Molecule 2
        # C atoms: [F,F,F,F,F,F,T,F,F,F,F]
        # Mol 2:   [F,F,F,F,F,F,T,T,T,T,T]
        # Expected: C6
        selector = (
            EnsembleSelector(ensemble, run_id).atom_types(["C"]).molecule_ids([2])
        )
        expected_mask = np.array(
            [
                False,
                False,
                False,
                False,
                False,
                False,
                True,
                False,
                False,
                False,
                False,
            ],
            dtype=np.bool_,
        )
        expected_masks = [expected_mask] * n_frames

        actual_masks_iterator = selector.get_mask(micro_id=slice(None))
        assert_masks_equal(actual_masks_iterator, expected_masks)

    def test_fluent_or_operation(
        self, test_project_with_synthetic_data: tuple[Project, Ensemble, str, int, int]
    ):
        project, ensemble, run_id, n_atoms, n_frames = test_project_with_synthetic_data

        # Select atoms that are (Hydrogen in Methane) OR (Oxygen)
        # H in Methane: [F,F,F,F,F,F,F,T,T,T,T]
        # Oxygen:       [T,F,F,T,F,F,F,F,F,F,F]
        # Expected: O0, O3, H7, H8, H9, H10
        selector = (
            EnsembleSelector(ensemble, run_id)
            .atom_types(
                ["H"]
            )  # This will pick up all 'H' and 'HW' (but in this case, only 'H' in methane)
            .OR()
            .atom_types(["OW"])
        )
        # Re-evaluating based on types in the fixture
        # atom_types = np.array(["OW", "HW", "HW", "OW", "HW", "HW", "C", "H", "H", "H", "H"])
        # atom_types(["H"]):                           [F,F,F,F,F,F,F,T,T,T,T]
        # atom_types(["OW"]):                          [T,F,F,T,F,F,F,F,F,F,F]
        # OR:                                            [T,F,F,T,F,F,F,T,T,T,T]
        expected_mask = np.array(
            [True, False, False, True, False, False, False, True, True, True, True],
            dtype=np.bool_,
        )
        expected_masks = [expected_mask] * n_frames

        actual_masks_iterator = selector.get_mask(micro_id=slice(None))
        assert_masks_equal(actual_masks_iterator, expected_masks)

    def test_fluent_not_operation(
        self, test_project_with_synthetic_data: tuple[Project, Ensemble, str, int, int]
    ):
        project, ensemble, run_id, n_atoms, n_frames = test_project_with_synthetic_data

        # Select atoms that are NOT "C" (Carbon)
        selector = EnsembleSelector(ensemble, run_id).NOT().atom_types(["C"])
        expected_mask = np.array(
            [True, True, True, True, True, True, False, True, True, True, True],
            dtype=np.bool_,
        )
        expected_masks = [expected_mask] * n_frames

        actual_masks_iterator = selector.get_mask(micro_id=slice(None))
        assert_masks_equal(actual_masks_iterator, expected_masks)

    def test_fluent_complex_chain(
        self, test_project_with_synthetic_data: tuple[Project, Ensemble, str, int, int]
    ):
        project, ensemble, run_id, n_atoms, n_frames = test_project_with_synthetic_data

        # Select (Oxygen OR Mol 2) AND (NOT Hydrogen)
        # O or Mol 2: [T,F,F,T,F,F,T,T,T,T,T]
        # NOT H:      [T,F,F,T,F,F,T,F,F,F,F] (original H values: [F,T,T,F,T,T,F,T,T,T,T])
        # AND result: [T,F,F,T,F,F,T,F,F,F,F] (O0, O3, C6)
        selector = (
            EnsembleSelector(ensemble, run_id)
            .atom_types(["OW"])
            .OR()
            .molecule_ids([2])
            .NOT()
            .atom_types(["HW", "H"])
        )
        expected_mask = np.array(
            [True, False, False, True, False, False, True, False, False, False, False],
            dtype=np.bool_,
        )
        expected_masks = [expected_mask] * n_frames

        actual_masks_iterator = selector.get_mask(micro_id=slice(None))
        assert_masks_equal(actual_masks_iterator, expected_masks)

    def test_fluent_implicit_and_vs_explicit_or(
        self, test_project_with_synthetic_data: tuple[Project, Ensemble, str, int, int]
    ):
        project, ensemble, run_id, n_atoms, n_frames = test_project_with_synthetic_data

        # Test implicit AND and explicit OR
        # (AtomTypeIs("OW") AND MolIdIs(0)) OR AtomTypeIs("C")
        # (O0, H1, H2) OR (C6) -> O0, H1, H2, C6
        selector = (
            EnsembleSelector(ensemble, run_id)
            .atom_types(["OW"])
            .molecule_ids([0])  # This implies an AND with the previous
            .OR()
            .atom_types(["C"])
        )
        # Breaking it down:
        # 1. ES().atom_types(["OW"]) -> O0, O3. Current: AtomTypeIs("OW")
        # 2. .molecule_ids([0]) -> Implicit AND. Current: And(AtomTypeIs("OW"), MolIdIs(0)) -> O0
        # 3. .OR() -> Sets next operation to OR
        # 4. .atom_types(["C"]) -> Current: Or(And(AtomTypeIs("OW"), MolIdIs(0)), AtomTypeIs("C"))
        # This simplifies to: (O0) OR (C6)
        expected_mask = np.array(
            [True, False, False, False, False, False, True, False, False, False, False],
            dtype=np.bool_,
        )
        expected_masks = [expected_mask] * n_frames

        actual_masks_iterator = selector.get_mask(micro_id=slice(None))
        assert_masks_equal(actual_masks_iterator, expected_masks)

    def test_get_mask_specific_microstate(
        self, test_project_with_synthetic_data: tuple[Project, Ensemble, str, int, int]
    ):
        project, ensemble, run_id, n_atoms, n_frames = test_project_with_synthetic_data

        # For microstate 0: Select atoms within 0.8 A of atom 0 (O0, (0,0,0))
        # Atom 0 itself is True. H1 (0.7,0.7,0) and H2 (0.7,-0.7,0) are both sqrt(0.98) ~ 0.99 away
        # So H1 and H2 are NOT within 0.8 A
        selector = EnsembleSelector(ensemble, run_id).distance_within(
            from_atoms=0, dist=0.8
        )
        expected_mask_ms0 = np.array(
            [
                True,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
            ],
            dtype=np.bool_,
        )

        actual_masks_iterator = selector.get_mask(micro_id=0)
        assert_masks_equal(actual_masks_iterator, [expected_mask_ms0])

    def test_get_mask_all_microstates_with_microstate_cadence_selector(
        self, test_project_with_synthetic_data: tuple[Project, Ensemble, str, int, int]
    ):
        project, ensemble, run_id, n_atoms, n_frames = test_project_with_synthetic_data

        # Select atoms within 1.0 A of atom 3 (O of second water)
        # Atom 3: (3.0, 0.0, 0.0) in MS0, (3.1, 0.0, 0.0) in MS1
        # H4: (3.7, 0.7, 0.0) in MS0, (3.8, 0.7, 0.0) in MS1
        # H5: (3.7, -0.7, 0.0) in MS0, (3.8, -0.7, 0.0) in MS1
        # Distances from Atom 3 (O3) to H4/H5 are sqrt(0.7^2+0.7^2) = 0.9899... < 1.0
        # So for both microstates, the mask is [F,F,F,T,T,T,F,F,F,F,F]
        selector = EnsembleSelector(ensemble, run_id).distance_within(
            from_atoms=3, dist=1.0
        )
        expected_mask = np.array(
            [False, False, False, True, True, True, False, False, False, False, False],
            dtype=np.bool_,
        )
        expected_masks = [expected_mask] * n_frames

        actual_masks_iterator = selector.get_mask(micro_id=slice(None))
        assert_masks_equal(actual_masks_iterator, expected_masks)

    def test_get_mask_missing_run_id(
        self,
        test_project_with_synthetic_data: tuple[Project, Ensemble, str, int, int],
        caplog,
    ):
        project, ensemble, run_id, n_atoms, n_frames = test_project_with_synthetic_data
        caplog.set_level(logging.WARNING)  # Expect warnings for missing data

        # Try to select from a non-existent run ID
        selector = EnsembleSelector(ensemble, "non_existent_run").atom_types(["OW"])

        # Mask should still work because atom_types are not on the run level.
        expected_mask = np.array(
            [True, False, False, True, False, False, False, False, False, False, False],
            dtype=np.bool_,
        )
        expected_masks = [
            expected_mask
        ] * n_frames  # Even if no frames, a warning should happen and then yield empty

        actual_masks_iterator = selector.get_mask(micro_id=slice(None))
        assert_masks_equal(actual_masks_iterator, expected_masks)

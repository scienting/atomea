# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

from typing import Any

import os
import random
import shutil

import numpy as np
import numpy.typing as npt
import polars as pl

from atomea.containers import Project
from atomea.stores import DiskFormat
from atomea.stores.arrays import ZarrArrayStore
from atomea.stores.tables import PolarsTableStore


class SyntheticAtomDataGenerator:
    """
    A robust and reusable class for generating synthetic atomistic data.

    This generator supports creating data for both deterministic testing (with fixed
    random seeds) and fuzzing (by injecting malformed data or using random seeds).
    """

    def __init__(
        self,
        n_frames: int | tuple[int, int] = (10, 50),
        n_atoms: int | tuple[int, int] = (50, 200),
        n_molecules: int | tuple[int, int] = (1, 10),
        n_components_per_mol: int | tuple[int, int] = (1, 5),
        random_seed: int | None = None,
        # Fuzzing control parameters: probabilities of injecting invalid data
        inject_nan_coords_prob: float = 0.0,
        inject_inf_coords_prob: float = 0.0,
        inject_invalid_shape_prob: float = 0.0,
        inject_invalid_dtype_prob: float = 0.0,
        inject_out_of_bounds_id_prob: float = 0.0,
        # Allow specific value ranges for fuzzing beyond typical values
        coord_min: float = -100.0,
        coord_max: float = 100.0,
    ):
        """
        Initializes the synthetic data generator.

        Args:
            n_frames: Number of frames. Can be a single int or a (min, max) tuple.
            n_atoms: Number of atoms. Can be a single int or a (min, max) tuple.
            n_molecules: Number of molecules. Can be a single int or a (min, max) tuple.
            n_components_per_mol: Number of components per molecule.
                Can be a single int or a (min, max) tuple.
            random_seed: If provided, seeds the random number generators for reproducible
                data generation. Set to None for truly random generation (fuzzing).
            inject_nan_coords_prob: Probability (0.0-1.0) of injecting NaN values into coordinates.
            inject_inf_coords_prob: Probability (0.0-1.0) of injecting Inf values into coordinates.
            inject_invalid_shape_prob: Probability (0.0-1.0) of generating arrays with
                incorrect shapes for fuzzing.
            inject_invalid_dtype_prob: Probability (0.0-1.0) of generating arrays/lists
                with incorrect data types for fuzzing.
            inject_out_of_bounds_id_prob: Probability (0.0-1.0) of injecting IDs
                that exceed typical limits (e.g., UInt32 max).
            coord_min: Minimum value for coordinates (for fuzzing, can be extreme).
            coord_max: Maximum value for coordinates (for fuzzing, can be extreme).
        """
        self.n_frames_range = (
            n_frames if isinstance(n_frames, tuple) else (n_frames, n_frames)
        )
        self.n_atoms_range = (
            n_atoms if isinstance(n_atoms, tuple) else (n_atoms, n_atoms)
        )
        self.n_molecules_range = (
            n_molecules
            if isinstance(n_molecules, tuple)
            else (n_molecules, n_molecules)
        )
        self.n_components_per_mol_range = (
            n_components_per_mol
            if isinstance(n_components_per_mol, tuple)
            else (n_components_per_mol, n_components_per_mol)
        )

        self.random_seed = random_seed
        if self.random_seed is not None:
            np.random.seed(self.random_seed)
            random.seed(self.random_seed)

        self.atom_symbols_pool = ["H", "C", "N", "O", "S", "P", "Fe", "Cl", "Na", "K"]
        self.atom_types_pool = ["HC", "CA", "OH", "N", "S", "FE", "CL", "NA", "K"]
        self.component_labels_pool = ["ALA", "GLY", "MET", "WAT", "ION", "LIG", "MEMB"]

        self.inject_nan_coords_prob = inject_nan_coords_prob
        self.inject_inf_coords_prob = inject_inf_coords_prob
        self.inject_invalid_shape_prob = inject_invalid_shape_prob
        self.inject_invalid_dtype_prob = inject_invalid_dtype_prob
        self.inject_out_of_bounds_id_prob = inject_out_of_bounds_id_prob
        self.coord_min = coord_min
        self.coord_max = coord_max

    def _get_random_size(self, size_range: tuple[int, int]) -> int:
        """Helper to get a random size within a specified range."""
        return random.randint(size_range[0], size_range[1])

    def generate_coordinates(
        self, n_frames: int, n_atoms: int
    ) -> npt.NDArray[np.float64]:
        """
        Generates synthetic atomic coordinates (n_frames, n_atoms, 3).
        Can inject NaNs, Infs, or malformed shapes/dtypes for fuzzing.
        """
        if random.random() < self.inject_invalid_shape_prob:
            # Inject incorrect shape for fuzzing
            if random.random() < 0.5:  # 2D array
                return np.random.rand(n_frames, n_atoms).astype(np.float64)
            else:  # 4D array
                return np.random.rand(n_frames, n_atoms, 3, 1).astype(np.float64)

        if random.random() < self.inject_invalid_dtype_prob:
            # Inject incorrect dtype for fuzzing (e.g., int instead of float)
            return np.random.randint(
                self.coord_min, self.coord_max, size=(n_frames, n_atoms, 3)
            ).astype(np.int32)

        coords = np.random.uniform(
            self.coord_min, self.coord_max, size=(n_frames, n_atoms, 3)
        ).astype(np.float64)

        if (
            self.inject_nan_coords_prob > 0
            and random.random() < self.inject_nan_coords_prob
        ):
            num_nans = max(1, int(n_frames * n_atoms * 3 * 0.01))  # Inject 1% NaNs
            for _ in range(num_nans):
                coords[
                    random.randint(0, n_frames - 1),
                    random.randint(0, n_atoms - 1),
                    random.randint(0, 2),
                ] = np.nan

        if (
            self.inject_inf_coords_prob > 0
            and random.random() < self.inject_inf_coords_prob
        ):
            num_infs = max(1, int(n_frames * n_atoms * 3 * 0.01))  # Inject 1% Infs
            for _ in range(num_infs):
                coords[
                    random.randint(0, n_frames - 1),
                    random.randint(0, n_atoms - 1),
                    random.randint(0, 2),
                ] = np.inf

        return coords

    def generate_atomic_numbers(self, n_atoms: int) -> npt.NDArray[np.uint8]:
        """
        Generates synthetic atomic numbers (n_atoms,).
        Can inject out-of-bounds values or malformed dtypes for fuzzing.
        """
        if random.random() < self.inject_invalid_dtype_prob:
            return np.random.uniform(1.0, 120.0, size=(n_atoms,)).astype(
                np.float32
            )  # Invalid dtype: float

        atomic_numbers = np.random.randint(
            1, 100, size=(n_atoms,), dtype=np.uint8
        )  # Max practical Z is 92 for U

        if (
            self.inject_out_of_bounds_id_prob > 0
            and random.random() < self.inject_out_of_bounds_id_prob
        ):
            if n_atoms > 0:
                # Inject a value > 255 for uint8
                idx = random.randint(0, n_atoms - 1)
                atomic_numbers[idx] = random.randint(
                    256, 500
                )  # This will cause overflow/wrap around for uint8
                # Alternative for true error if not using numpy's strict type enforcement:
                # return atomic_numbers.astype(np.int16) # Return a dtype that doesn't fit uint8

        return atomic_numbers

    def generate_atom_symbols(self, n_atoms: int) -> list[str]:
        """
        Generates synthetic elemental symbols (list of strings).
        Can inject invalid types or incorrect counts for fuzzing.
        """
        if random.random() < self.inject_invalid_dtype_prob:
            return [
                random.randint(1, 100) for _ in range(n_atoms)
            ]  # Invalid type: int list

        if random.random() < self.inject_invalid_shape_prob:
            # Return list with incorrect number of elements
            return random.choices(
                self.atom_symbols_pool, k=n_atoms + random.randint(-5, 5)
            )

        return random.choices(self.atom_symbols_pool, k=n_atoms)

    def generate_atom_types(self, n_atoms: int) -> list[str]:
        """
        Generates synthetic classical force field atom types (list of strings).
        """
        # Minimal fuzzing for now, can extend with invalid/malformed strings later
        return random.choices(self.atom_types_pool, k=n_atoms)

    def generate_molecule_ids(
        self, n_atoms: int, n_molecules: int
    ) -> npt.NDArray[np.uint32]:
        """
        Generates synthetic molecule IDs (n_atoms,).
        Can inject out-of-bounds values or malformed dtypes/shapes.
        """
        if random.random() < self.inject_invalid_shape_prob:
            return np.random.randint(0, n_molecules, size=(n_atoms, 2)).astype(
                np.uint32
            )  # 2D array

        if random.random() < self.inject_invalid_dtype_prob:
            return np.random.uniform(0.0, n_molecules, size=(n_atoms,)).astype(
                np.float64
            )  # Invalid dtype: float

        ids = np.zeros(n_atoms, dtype=np.uint32)
        if n_molecules > 0 and n_atoms > 0:
            atoms_per_mol = max(1, n_atoms // n_molecules)
            current_mol_id = 0
            for i in range(n_atoms):
                ids[i] = current_mol_id
                if (i + 1) % atoms_per_mol == 0 and current_mol_id < n_molecules - 1:
                    current_mol_id += 1

        if (
            self.inject_out_of_bounds_id_prob > 0
            and random.random() < self.inject_out_of_bounds_id_prob
        ):
            if n_atoms > 0:
                # Inject an ID too large for UInt32
                idx = random.randint(0, n_atoms - 1)
                ids[idx] = 2**32 + random.randint(
                    1, 1000
                )  # This will cause overflow/error for UInt32

        return ids

    def generate_component_ids(
        self, n_atoms: int, n_components_per_mol: int
    ) -> npt.NDArray[np.uint32]:
        """
        Generates synthetic component IDs (n_atoms,).
        Can inject out-of-bounds values or malformed dtypes/shapes.
        """
        if random.random() < self.inject_invalid_shape_prob:
            return (
                np.random.randint(0, n_components_per_mol, size=(n_atoms,))
                .reshape(-1, 1)
                .astype(np.uint32)
            )  # Column vector

        if random.random() < self.inject_invalid_dtype_prob:
            return np.random.uniform(0.0, n_components_per_mol, size=(n_atoms,)).astype(
                np.int64
            )  # Invalid dtype: int64

        ids = np.zeros(n_atoms, dtype=np.uint32)
        if n_atoms > 0:
            current_comp_id = 1  # Start from 1 as per your example
            for i in range(n_atoms):
                ids[i] = current_comp_id
                # Increment component ID after a few atoms, or based on a probability
                if random.random() < (1 / max(1, n_components_per_mol)):
                    current_comp_id += 1

        if (
            self.inject_out_of_bounds_id_prob > 0
            and random.random() < self.inject_out_of_bounds_id_prob
        ):
            if n_atoms > 0:
                idx = random.randint(0, n_atoms - 1)
                ids[idx] = 2**32 + random.randint(1, 1000)

        return ids

    def generate_component_labels(self, n_atoms: int) -> list[str]:
        """
        Generates synthetic component labels (list of strings).
        """
        # Minimal fuzzing for now, can extend with invalid/malformed strings later
        return random.choices(self.component_labels_pool, k=n_atoms)

    def generate_energy_data(
        self, n_frames: int, ens_id: str, run_id: str
    ) -> pl.DataFrame:
        """
        Generates synthetic energy data as a Polars DataFrame.
        This includes 'ens_id', 'run_id', 'micro_id', and 'potential_energy'.
        """
        # For fuzzing, we could inject missing columns, wrong dtypes in columns, etc.
        # For now, generate valid data as Polars handles many data issues.
        data = {
            "ens_id": [ens_id] * n_frames,
            "run_id": [run_id] * n_frames,
            "micro_id": list(range(n_frames)),
            "potential_mm": np.random.uniform(-500.0, 500.0, size=n_frames).astype(
                np.float64
            ),
        }
        return pl.DataFrame(data)

    def create_synthetic_project(
        self,
        tmp_dir: str,
        ens_id: str = "synthetic_ensemble",
        run_id: str = "run_0",
        clean_up_paths: bool = True,
    ) -> tuple[Project, dict[str, Any]]:
        """
        Generates synthetic data and writes it to a new Project instance.

        Args:
            tmp_dir: Temporary directory path where Zarr and Polars stores will be created.
            ens_id: ID for the ensemble to be created.
            run_id: ID for the run within the ensemble.
            clean_up_paths: If True, remove existing store directories before creating.

        Returns:
            A tuple containing:
            - The initialized Project instance with synthetic data written.
            - A dictionary of the generated reference data (before writing), useful for assertions.
        """
        n_frames = self._get_random_size(self.n_frames_range)
        n_atoms = self._get_random_size(self.n_atoms_range)
        n_molecules = self._get_random_size(self.n_molecules_range)
        n_components_per_mol = self._get_random_size(self.n_components_per_mol_range)

        # Generate all data first
        ref_coords = self.generate_coordinates(n_frames, n_atoms)
        ref_atomic_numbers = self.generate_atomic_numbers(n_atoms)
        ref_symbols = self.generate_atom_symbols(n_atoms)
        ref_types = self.generate_atom_types(n_atoms)
        ref_mol_ids = self.generate_molecule_ids(n_atoms, n_molecules)
        ref_comp_ids = self.generate_component_ids(n_atoms, n_components_per_mol)
        ref_comp_labels = self.generate_component_labels(n_atoms)
        ref_energy_df = self.generate_energy_data(n_frames, ens_id, run_id)

        path_store_array = os.path.join(tmp_dir, f"{ens_id}_{run_id}_arrays.zarr")
        path_store_table = os.path.join(tmp_dir, f"{ens_id}_{run_id}_tables.tables")

        if clean_up_paths:
            if os.path.exists(path_store_array):
                shutil.rmtree(path_store_array)
            if os.path.exists(path_store_table):
                shutil.rmtree(path_store_table)

        # Initialize stores and project
        store_array = ZarrArrayStore(path_store_array, mode="a")
        store_table = PolarsTableStore(
            path_store_table, mode="a", disk_format=DiskFormat.PARQUET
        )
        prj = Project(store_array, store_table)
        ensemble = prj.add_ensemble(ens_id)

        # Write data to the project. This is where validation errors are expected during fuzzing.
        ensemble.coordinates.write(ref_coords, run_id=run_id)
        ensemble.topology.atoms.atomic_numbers.write(ref_atomic_numbers, run_id=run_id)
        ensemble.topology.atoms.symbols.write(np.array(ref_symbols), run_id=run_id)
        ensemble.topology.atoms.types.write(np.array(ref_types), run_id=run_id)
        ensemble.topology.ids.molecules.write(ref_mol_ids, run_id=run_id)
        ensemble.topology.ids.components.write(ref_comp_ids, run_id=run_id)
        ensemble.topology.labels.components.write(
            np.array(ref_comp_labels), run_id=run_id
        )
        prj.energy.potential_mm.write(ref_energy_df, run_id=run_id)

        reference_data = {
            "n_frames": n_frames,
            "n_atoms": n_atoms,
            "ref_coords": ref_coords,
            "ref_atomic_numbers": ref_atomic_numbers,
            "ref_symbols": ref_symbols,
            "ref_types": ref_types,
            "ref_mol_ids": ref_mol_ids,
            "ref_comp_ids": ref_comp_ids,
            "ref_comp_labels": ref_comp_labels,
            "ref_energy_df": ref_energy_df,
            "ens_id": ens_id,
            "run_id": run_id,
        }
        return prj, reference_data

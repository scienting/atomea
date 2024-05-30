"""Processing AMBER simulations with roGFP2"""

import os

import numpy as np

from atomea.digesters import MDAnalysisDigester

from .conftest import TMP_DIR

# from atomea.io.backends._zarr import ZarrStorage
# from atomea.io.write import DataProcessor


def test_digest_amber_rogfp2_serial(amber_rogfp2_sim_paths):
    """
    Test the digestion of Amber simulations using the roGFP2 ensemble.

    This test function verifies the correct digestion of Amber simulations
    using the roGFP2 ensemble. It uses the `EnsembleSchema` class to store
    the ensemble data and the `MDAnalysisDigester` class to perform the
    digestion. The digestion is performed by passing the necessary input
    files, such as the topology file and the trajectory file, to the digester.
    """
    digester = MDAnalysisDigester()
    digester_args = (
        amber_rogfp2_sim_paths["mol.prmtop"],
        amber_rogfp2_sim_paths["07_relax_npt.nc"],
    )
    digester_kwargs = {"topology_format": "PRMTOP", "format": "NC"}
    ensemble_data = digester.digest(
        digester_args=digester_args, digester_kwargs=digester_kwargs
    )
    assert len(ensemble_data.frames) == 100
    assert ensemble_data.frames[0].system.coordinates is not None
    assert np.allclose(ensemble_data.frames[0].system.coordinates[0][0], 33.924496)
    assert np.allclose(ensemble_data.frames[0].system.coordinates[32][0], 27.2496)
    assert ensemble_data.frames[-1].system.coordinates is not None
    assert np.allclose(ensemble_data.frames[-1].system.coordinates[78][0], 29.406982)


# def test_write_amber_rogfp2(amber_rogfp2_sim_paths):
#     data_processor = DataProcessor(
#         digester_cls=MDAnalysisDigester,
#         storage_cls=ZarrStorage,
#         chunk_size=100,
#     )
#     zarr_path = os.path.join(TMP_DIR, "amber_rogfp2.zarr")
#     if os.path.exists(zarr_path):
#         os.remove(zarr_path)

#     digester_args = (
#         amber_rogfp2_sim_paths["mol.prmtop"],
#         amber_rogfp2_sim_paths["07_relax_npt.nc"],
#     )
#     digester_kwargs = {"topology_format": "PRMTOP", "format": "NC"}

#     data_processor.process_data(
#         zarr_path,
#         digester_args=digester_args,
#         digester_kwargs=digester_kwargs,
#         parallelism=4,
#     )

#     assert os.path.exists(zarr_path)

"""Processing AMBER simulations with roGFP2"""
from atomea.digesters import MDAnalysisDigester
from atomea.schemas import EnsembleSchema


def test_amber_rogfp2_digest(amber_rogfp2_sim_paths):
    """
    Test the digestion of Amber simulations using the roGFP2 ensemble.

    This test function verifies the correct digestion of Amber simulations
    using the roGFP2 ensemble. It uses the `EnsembleSchema` class to store
    the ensemble data and the `MDAnalysisDigester` class to perform the
    digestion. The digestion is performed by passing the necessary input
    files, such as the topology file and the trajectory file, to the digester.
    """
    ensemble_data = EnsembleSchema()
    digester = MDAnalysisDigester()
    digest_inputs = {"topology_format": "PRMTOP", "format": "NC"}
    ensemble_data = digester.digest(
        ensemble_data,
        amber_rogfp2_sim_paths["mol.prmtop"],
        amber_rogfp2_sim_paths["07_relax_npt.nc"],
        **digest_inputs
    )

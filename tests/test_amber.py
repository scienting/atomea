"""Processing AMBER simulations with roGFP2"""
from atomea.digesters import MDAnalysisDigester
from atomea.schemas import EnsembleSchema


def test_amber_rogfp2_digest(amber_rogfp2_sim_paths):
    ensemble_data = EnsembleSchema()
    digester = MDAnalysisDigester()
    digest_inputs = {
        "topology": amber_rogfp2_sim_paths["mol.prmtop"],
        "coordinates": amber_rogfp2_sim_paths["07_relax_npt.nc"],
    }
    ensemble_data = digester.digest(ensemble_data, **digest_inputs)

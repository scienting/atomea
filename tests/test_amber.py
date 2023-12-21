"""Processing AMBER simulations with roGFP2"""
import os
import tempfile

from atomea.digesters import MDAnalysisDigester
from atomea.io import DiskData


def test_amber_rogfp2_digest(uuid_simlify_rogfp2, base_atomea):
    digester = MDAnalysisDigester()
    topo_path = os.path.join(uuid_simlify_rogfp2, "topo/mol.prmtop")
    coord_path = os.path.join(uuid_simlify_rogfp2, "outputs/07_relax_npt.nc")

    schema = base_atomea.get()
    keep_keys = "coordinates"
    base_atomea.schema = {key: schema[key] for key in schema if key in keep_keys}
    data = digester.digest(base_atomea, topo_path, coord_path)

    disk = DiskData()
    with tempfile.TemporaryDirectory() as tmp_dir:
        disk.store(tmp_dir, data, base_atomea)

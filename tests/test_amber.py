"""Processing AMBER simulations with roGFP2"""
import os
import tempfile

import pyarrow as pa
import pyarrow.parquet as pq

from atomea.digesters import MDAnalysisDigester
from atomea.io import DiskData


def test_amber_rogfp2_digest(uuid_simlify_rogfp2, base_atomea):
    digester = MDAnalysisDigester()
    topo_path = os.path.join(uuid_simlify_rogfp2, "topo/mol.prmtop")
    coord_path = os.path.join(uuid_simlify_rogfp2, "outputs/07_relax_npt.nc")

    base_atomea.filter(("coordinates", "ff_atom_type"))
    data = digester.digest(base_atomea, topo_path, coord_path)

    disk = DiskData()
    with tempfile.TemporaryDirectory() as tmp_dir:
        disk.store(tmp_dir, data, base_atomea)

        table_test = pq.read_table(os.path.join(tmp_dir, "atoms.parquet"))

        assert table_test.column("ff_atom_type").type.equals(pa.string())
        assert table_test.shape == (373810, 1)

import os
import shutil

from atomea.containers import Project
from atomea.digesters.amber import AmberOutputDigester
from atomea.digesters.amber.v22 import (
    AmberV22Parser,
)
from atomea.stores import DiskFormat
from atomea.stores.arrays import ZarrArrayStore
from atomea.stores.tables import PolarsTableStore


class TestOutputDigestion:
    """Test digestion of amber output files"""

    def test_amber_npt_output(self, amber_rogfp2_sim_paths, tmp_dir):
        path_store_array = os.path.join(tmp_dir, "amber_rogfp2_npt_output.zarr")
        if os.path.exists(path_store_array):
            shutil.rmtree(path_store_array)

        path_store_table = os.path.join(tmp_dir, "amber_rogfp2_npt_output.tables")
        if os.path.exists(path_store_table):
            shutil.rmtree(path_store_table)

        store_array = ZarrArrayStore(path_store_array, mode="a")
        store_table = PolarsTableStore(path_store_table, DiskFormat.CSV)

        prj = Project(store_array, store_table)
        parser = AmberV22Parser()
        digest_args = (
            amber_rogfp2_sim_paths["07_relax_npt.out"],
            parser,
        )
        digest_kwargs = {}
        prj = AmberOutputDigester.run(
            prj,
            "default",
            "0",
            digest_args,
            digest_kwargs,
        )
        # Duplicate to test appending
        prj = AmberOutputDigester.run(
            prj,
            "default",
            "0",
            digest_args,
            digest_kwargs,
        )

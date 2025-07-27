# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

import os
import shutil

from atomea.containers import Project
from atomea.io.formats.amber import AmberOutputReader
from atomea.io.formats.amber.v22 import (
    AmberV22Parser,
)
from atomea.stores import DiskFormat
from atomea.stores.arrays import ZarrArrayStore
from atomea.stores.tables import PolarsTableStore


class TestOutputReading:
    """Test reading of amber output files"""

    def test_amber_npt_output(self, amber_rogfp2_sim_paths, tmp_dir):
        path_store_array = os.path.join(tmp_dir, "amber_rogfp2_npt_output.zarr")
        if os.path.exists(path_store_array):
            shutil.rmtree(path_store_array)

        path_store_table = os.path.join(tmp_dir, "amber_rogfp2_npt_output.tables")
        if os.path.exists(path_store_table):
            shutil.rmtree(path_store_table)

        store_array = ZarrArrayStore(path_store_array, mode="a")
        store_table = PolarsTableStore(
            path_store_table, mode="a", disk_format=DiskFormat.CSV
        )

        prj = Project(store_array, store_table)
        parser = AmberV22Parser()
        reader_args = (
            amber_rogfp2_sim_paths["07_relax_npt.out"],
            parser,
        )
        reader_kwargs = {}
        prj = AmberOutputReader.run(
            prj,
            "default",
            "0",
            reader_args,
            reader_kwargs,
        )
        # Duplicate to test appending
        prj = AmberOutputReader.run(
            prj,
            "default",
            "0",
            reader_args,
            reader_kwargs,
        )

import os
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import List


@dataclass
class SetupBinpath:
    """Setup binary path Class"""

    root_binpath: Path
    bin_list: List[str] = field(
        default_factory=lambda: [
            "make_ultrametric.py",
            "orthofinder.py",
            "mafft",
            "trimal",
            "iqtree",
            "AnGST.py",
            "AnGST_wrapper.py",
            "parallel",
        ]
    )

    def __post_init__(self):
        self._add_bin_path()
        self._bin_exists_check(self.bin_list)

    def _add_bin_path(self) -> None:
        """Add bin programs path"""
        # Bin path list
        mafft_path = self.root_binpath / "mafft"
        ortho_finder_path = self.root_binpath / "OrthoFinder"
        ortho_finder_tool_path = ortho_finder_path / "tools"
        angst_path = self.root_binpath / "angst" / "angst_lib"
        # Add bin path
        env_path = os.environ["PATH"]
        env_path = f"{self.root_binpath}:{env_path}"
        env_path = f"{mafft_path}:{env_path}"
        env_path = f"{ortho_finder_path}:{env_path}"
        env_path = f"{ortho_finder_tool_path}:{env_path}"
        env_path = f"{angst_path}:{env_path}"
        # Set fixed path
        os.environ["PATH"] = env_path

    def _bin_exists_check(self, bin_list: List[str]) -> None:
        """Check bin program exists"""
        bin_exists_check_flg = False
        print("Required bin program exists check...")
        for bin in bin_list:
            bin_path = shutil.which(bin)
            if bin_path:
                print(f"OK '{bin}' (Path={bin_path})")
            else:
                print(f"NG '{bin}' (Path=Not exists)")
                bin_exists_check_flg = True
        if bin_exists_check_flg:
            print("Required bin program not exist!!\n")
            exit(1)
        else:
            print("All required bin programs exists.\n")

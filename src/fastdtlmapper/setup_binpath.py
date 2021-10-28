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
            "treerecs",
            "AnGST.py",
            "AnGST_wrapper.py",
            "parallel",
        ]
    )

    def __post_init__(self):
        self._add_bin_path(self._get_bin_path_list())
        self._bin_exists_check(self.bin_list)

    def _get_bin_path_list(self) -> List[Path]:
        return [
            self.root_binpath,
            self.root_binpath / "mafft",
            self.root_binpath / "OrthoFinder",
            self.root_binpath / "OrthoFinder" / "tools",
            self.root_binpath / "angst" / "angst_lib",
        ]

    def _add_bin_path(self, bin_path_list: List[Path]) -> None:
        """Add bin programs path"""
        # Add bin path
        env_path = os.environ["PATH"]
        for bin_path in bin_path_list:
            env_path = f"{bin_path}:{env_path}"
        # Set fixed path
        os.environ["PATH"] = env_path

    def _bin_exists_check(self, bin_list: List[str]) -> None:
        """Check bin program exists"""
        bin_exists_check_flg = False
        print("Checking required bin programs...")
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
            print("Check result = All OK!!\n")

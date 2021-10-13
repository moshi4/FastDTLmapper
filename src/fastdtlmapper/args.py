from __future__ import annotations

import argparse
import datetime as dt
import os
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Union


@dataclass
class Args:
    """Parse Arguments Class"""

    indir: Path
    tree_file: Path
    outdir: Path
    process_num: int
    dup_cost: int
    los_cost: int
    trn_cost: int
    inflation: float
    timetree: bool
    rseed: int

    def __post_init__(self):
        self._start_time = time.time()

    def write_log(self, log_file: Union[str, Path]) -> None:
        """Write arguments log file

        Args:
            log_file (Union[str, Path]): Output log file path
        """
        with open(log_file, "w") as f:
            f.write(self.log_text)

    @property
    def log_text(self) -> str:
        """Get arguments log"""
        end_time = time.time()
        elapsed_time = (end_time - self._start_time) / 3600

        def unixtime_to_datestr(unixtime: float) -> str:
            return dt.datetime.fromtimestamp(unixtime).strftime("%Y/%m/%y %H:%M:%S")

        log_text = (
            "Run command:\n"
            + f"{' '.join(sys.argv)}\n\n"
            + f"Input directory = {self.indir}\n"
            + f"Tree file = {self.tree_file}\n"
            + f"Output directory = {self.outdir}\n"
            + f"Number of processor = {self.process_num}\n"
            + f"Duplication event cost = {self.dup_cost}\n"
            + f"Loss event cost = {self.los_cost}\n"
            + f"Transfer event cost = {self.trn_cost}\n"
            + f"MCL inflation parameter = {self.inflation}\n"
            + f"Number of random seed = {self.rseed}\n\n"
            + f"Start time = {unixtime_to_datestr(self._start_time)}\n"
            + f"End time = {unixtime_to_datestr(end_time)}\n"
            + f"Elapsed time = {elapsed_time:.2f}[h]"
        )
        return log_text


def get_args(argv: Optional[List[str]] = None) -> Args:
    """Get arguments

    Returns:
        Args: Args Class
    """
    parser = argparse.ArgumentParser(
        description="Fast genome-wide DTL event mapping tool"
    )

    parser.add_argument(
        "-i",
        "--indir",
        required=True,
        type=Path,
        help="Input Fasta(*.fa|*.faa|*.fasta), "
        + "Genbank(*.gb|*.gbk|*.genbank) directory",
        metavar="IN",
    )
    parser.add_argument(
        "-t",
        "--tree",
        required=True,
        type=Path,
        help="Input rooted species newick tree file (timetree is preferable)",
        metavar="TREE",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        required=True,
        type=Path,
        help="Output directory",
        metavar="OUT",
    )
    cpu_count = os.cpu_count()
    default_processor_num = cpu_count - 1 if cpu_count is not None else 1
    parser.add_argument(
        "-p",
        "--process_num",
        type=int,
        help=f"Number of processor (Default: {default_processor_num})",
        default=default_processor_num,
        metavar="",
    )
    default_dup_cost = 2
    parser.add_argument(
        "--dup_cost",
        type=int,
        help=f"Duplication event cost (Default: {default_dup_cost})",
        default=default_dup_cost,
        metavar="",
    )
    default_los_cost = 1
    parser.add_argument(
        "--los_cost",
        type=int,
        help=f"Loss event cost (Default: {default_los_cost})",
        default=default_los_cost,
        metavar="",
    )
    default_trn_cost = 3
    parser.add_argument(
        "--trn_cost",
        type=int,
        help=f"Transfer event cost (Default: {default_trn_cost})",
        default=default_trn_cost,
        metavar="",
    )
    default_inflation = 3.0
    parser.add_argument(
        "--inflation",
        type=float,
        help=f"MCL inflation parameter (Default: {default_inflation})",
        default=default_inflation,
        metavar="",
    )
    parser.add_argument(
        "--timetree",
        help="Use species tree as timetree",
        action="store_true",
    )
    default_rseed = 0
    parser.add_argument(
        "--rseed",
        type=int,
        help=f"Number of random seed (Default: {default_rseed})",
        default=default_rseed,
        metavar="",
    )

    args = parser.parse_args(argv)

    return Args(
        args.indir,
        args.tree,
        args.outdir,
        args.process_num,
        args.dup_cost,
        args.los_cost,
        args.trn_cost,
        args.inflation,
        args.timetree,
        args.rseed,
    )

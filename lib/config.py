import argparse
import os
import random
import shutil
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import List


@dataclass
class Config:
    """Config Class"""

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
        self.start_time = time.time()
        # 00. User fasta & tree directory
        self.user_data_dir = self.outdir / "00_user_data"
        self.user_fasta_dir = self.user_data_dir / "fasta"
        self.user_tree_dir = self.user_data_dir / "tree"

        self.user_tree_file = self.user_tree_dir / "user_tree.nwk"
        self.um_tree_file = self.user_tree_dir / "ultrametric_tree.nwk"
        self.um_nodeid_tree_file = self.user_tree_dir / "ultrametric_nodeid_tree.nwk"

        # 01. OrthoFinder directory
        self.ortho_dir = self.outdir / "01_orthofinder"
        self.ortho_group_fasta_dir = self.ortho_dir / "Orthogroup_Sequences"
        self.ortho_tmpwork_dir = self.user_fasta_dir / "OrthoFinder"

        # 02. DTL reconciliation directory
        self.dtl_rec_dir = self.outdir / "02_dtl_reconciliation"

        # 03. All group node DTL event aggregate map directory
        self.aggregate_map_dir = self.outdir / "03_aggregate_map_result"

        self.all_dtl_map_file = self.aggregate_map_dir / "all_dtl_map.nwk"
        self.all_gain_loss_map_file = self.aggregate_map_dir / "all_gain_loss_map.nwk"
        self.all_node_event_file = self.aggregate_map_dir / "all_group_node_event.tsv"
        self.all_transfer_gene_file = self.aggregate_map_dir / "all_transfer_genes.tsv"

        # Log commands list and command stderr
        self.log_dir = self.outdir / "log"
        self.parallel_cmds_log_dir = self.log_dir / "parallel_cmds"
        self.cmd_stderr_log_dir = self.log_dir / "cmds_stderr"

        self.run_config_log_file = self.log_dir / "run_config.log"

        self.mafft_stderr_log_file = self.cmd_stderr_log_dir / "mafft_stderr.log"
        self.trimal_stderr_log_file = self.cmd_stderr_log_dir / "trimal_stderr.log"
        self.iqtree_stderr_log_file = self.cmd_stderr_log_dir / "iqtree_stderr.log"

        self._makedirs()
        self._delete_prev_stderr_log()
        self._add_bin_path()
        self._bin_exists_check()

    def _makedirs(self) -> None:
        """Create config output directories"""
        # Check & Delete OrthoFinder previous work dir
        if self.ortho_tmpwork_dir.exists():
            shutil.rmtree(self.ortho_tmpwork_dir)

        # Make configured output directory
        for target_dir in (
            self.user_fasta_dir,
            self.user_tree_dir,
            self.ortho_dir,
            self.dtl_rec_dir,
            self.aggregate_map_dir,
            self.parallel_cmds_log_dir,
            self.cmd_stderr_log_dir,
        ):
            os.makedirs(target_dir, exist_ok=True)

    def _delete_prev_stderr_log(self) -> None:
        """Delete previous stderr log"""
        for stderr_log_file in (
            self.mafft_stderr_log_file,
            self.trimal_stderr_log_file,
            self.iqtree_stderr_log_file,
        ):
            stderr_log_file.unlink(missing_ok=True)

    def _add_bin_path(self) -> None:
        """Add bin programs path"""
        # Bin path list
        bin_path = Path(__file__).parent.parent / "bin"
        mafft_path = bin_path / "mafft"
        ortho_finder_path = bin_path / "OrthoFinder"
        ortho_finder_tool_path = ortho_finder_path / "tools"
        angst_path = bin_path / "angst" / "angst_lib"
        # Add bin path
        env_path = os.environ["PATH"]
        env_path = f"{bin_path}:{env_path}"
        env_path = f"{mafft_path}:{env_path}"
        env_path = f"{ortho_finder_path}:{env_path}"
        env_path = f"{ortho_finder_tool_path}:{env_path}"
        env_path = f"{angst_path}:{env_path}"
        # Set fixed path
        os.environ["PATH"] = env_path

    def _bin_exists_check(self) -> None:
        """Check bin program exists"""
        bin_list = [
            "make_ultrametric.py",
            "orthofinder.py",
            "mafft",
            "trimal",
            "iqtree",
            "AnGST.py",
            "AnGST_wrapper.py",
            "parallel",
        ]
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

    def output_run_config_log(self):
        """Output run config log file"""
        run_command = " ".join(sys.argv)
        elapsed_time = (time.time() - self.start_time) / 3600

        output_log = ""
        output_log += "Run command:\n"
        output_log += f"{run_command}\n"
        output_log += "\n"
        output_log += f"Input directory = {self.indir}\n"
        output_log += f"Tree file = {self.tree_file}\n"
        output_log += f"Output directory = {self.outdir}\n"
        output_log += f"Number of processor = {self.process_num}\n"
        output_log += f"Duplication event cost = {self.dup_cost}\n"
        output_log += f"Loss event cost = {self.los_cost}\n"
        output_log += f"Transfer event cost = {self.trn_cost}\n"
        output_log += f"MCL inflation parameter = {self.inflation}\n"
        output_log += f"Number of random seed = {self.rseed}\n"
        output_log += "\n"
        output_log += f"Elapsed time = {elapsed_time:.2f}[h]"

        with open(self.run_config_log_file, "w") as f:
            f.write(output_log)

    def make_ultrametric_cmd(self, nwk_tree_infile: str, nwk_tree_outfile: str) -> str:
        """Get make_ultrametric run command"""
        return f"make_ultrametric.py {nwk_tree_infile} {nwk_tree_outfile} "

    def orthofinder_cmd(self, fasta_indir: str) -> str:
        """Get OrthoFinder run command"""
        return (
            f"orthofinder.py -og -f {fasta_indir} -t {self.process_num} "
            + f"-I {self.inflation}"
        )

    def mafft_cmd(self, fasta_infile: str, aln_outfile: str) -> str:
        """Get mafft run command"""
        return (
            f"mafft --auto --anysymbol --quiet {fasta_infile} > {aln_outfile} "
            + f"2>> {self.mafft_stderr_log_file}"
        )

    def trimal_cmd(self, aln_infile: str, aln_trim_outfile: str) -> str:
        """Get trimal run command"""
        return (
            f"trimal -in {aln_infile} -out {aln_trim_outfile} -automated1 "
            + f"2>> {self.trimal_stderr_log_file}"
        )

    def iqtree_cmd(self, aln_infile: str, prefix: str) -> str:
        """Get iqtree run command"""
        # TODO: Set model parameter from argument
        return (
            f"iqtree -s {aln_infile} --prefix {prefix} -m TEST -mset JTT,WAG,LG "
            + f"--seed {self.rseed} --ufboot 1000 --boot-trees --wbtl --redo --quiet "
            + f"2>> {self.iqtree_stderr_log_file} "
        )

    def angst_cmd(self, species_tree_file: str, boot_tree_file: str, outdir: str):
        """Get AnGST run command"""
        timetree = "--timetree" if self.timetree else ""
        return (
            f"AnGST_wrapper.py -s {species_tree_file} -g {boot_tree_file} -o {outdir} "
            + f"--dup_cost {self.dup_cost} --los_cost {self.los_cost} "
            + f"--trn_cost {self.trn_cost} {timetree} "
        )

    def parallel_cmd(self, cmd_list: List[str], prog_name: str) -> None:
        """Get parallel run command from command list"""
        random.seed(self.rseed)
        random.shuffle(cmd_list)
        parallel_cmds_file = self.parallel_cmds_log_dir / (prog_name + "_cmds.log")
        with open(parallel_cmds_file, "w") as f:
            for cmd in cmd_list:
                f.write(cmd + "\n")
        return (
            f"parallel --no-notice --bar -a {parallel_cmds_file} "
            + f"-j {self.process_num} "
        )


def get_config() -> Config:
    """Get config from arguments

    Returns:
        Config: Config Class
    """
    parser = argparse.ArgumentParser(
        description="Fast genome-wide DTL event mapping tool"
    )

    parser.add_argument(
        "-i",
        "--indir",
        required=True,
        type=Path,
        help="Input Fasta(*.fa|*.faa|*.fasta), Genbank(*.gb|*.gbk|*.genbank) directory",
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
    default_processor_num = os.cpu_count() - 1
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

    args = parser.parse_args()

    return Config(
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

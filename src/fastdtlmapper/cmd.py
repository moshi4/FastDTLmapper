import random
import subprocess as sp
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Union

from fastdtlmapper.args import Args, get_args


@dataclass
class Cmd:
    """Run Command Class"""

    params: Args = field(default_factory=get_args)

    def get_make_ultrametric_cmd(
        self, nwk_tree_infile: Union[str, Path], nwk_tree_outfile: Union[str, Path]
    ) -> str:
        """Get make_ultrametric run command"""
        return f"make_ultrametric.py {nwk_tree_infile} {nwk_tree_outfile} "

    def get_orthofinder_cmd(self, fasta_indir: Union[str, Path]) -> str:
        """Get OrthoFinder run command"""
        return (
            f"orthofinder.py -og -f {fasta_indir} -t {self.params.process_num} "
            + f"-I {self.params.inflation}"
        )

    def get_mafft_cmd(
        self, fasta_infile: Union[str, Path], aln_outfile: Union[str, Path]
    ) -> str:
        """Get mafft run command"""
        return f"mafft --auto --anysymbol --quiet {fasta_infile} > {aln_outfile} 2>&1"

    def get_trimal_cmd(
        self, aln_infile: Union[str, Path], aln_trim_outfile: Union[str, Path]
    ) -> str:
        """Get trimal run command"""
        return f"trimal -in {aln_infile} -out {aln_trim_outfile} -automated1 2>&1"

    def get_iqtree_cmd(
        self, aln_infile: Union[str, Path], prefix: Union[str, Path], boot: bool = True
    ) -> str:
        """Get iqtree run command"""
        boot_opt = "--ufboot 1000 --boot-trees --wbtl" if boot else ""
        return (
            f"iqtree -s {aln_infile} --prefix {prefix} -m TEST -mset JTT,WAG,LG "
            + f"--seed {self.params.rseed} {boot_opt} --redo --quiet 2>&1"
        )

    def get_angst_cmd(
        self,
        species_tree_file: Union[str, Path],
        boot_tree_file: Union[str, Path],
        outdir: Union[str, Path],
    ) -> str:
        """Get AnGST run command"""
        timetree = "--timetree" if self.params.timetree else ""
        return (
            f"AnGST_wrapper.py -s {species_tree_file} -g {boot_tree_file} -o {outdir} "
            + f"--dup_cost {self.params.dup_cost} --los_cost {self.params.los_cost} "
            + f"--trn_cost {self.params.trn_cost} {timetree} 2>&1"
        )

    def get_parallel_cmd(
        self,
        cmd_list: List[str],
        parallel_cmds_file: Union[str, Path],
        parallel_log_file: Union[str, Path],
    ) -> str:
        """Get parallel run command from command list"""
        random.seed(self.params.rseed)
        random.shuffle(cmd_list)
        with open(parallel_cmds_file, "w") as f:
            for cmd in cmd_list:
                f.write(cmd + "\n")
        return (
            f"parallel --no-notice --bar -a {parallel_cmds_file} "
            + f"-j {self.params.process_num} --results {parallel_log_file} "
            + "> /dev/null "
        )

    def run_parallel_cmd(
        self,
        cmd_list: List[str],
        parallel_cmds_file: Union[str, Path],
        parallel_log_file: Union[str, Path],
    ) -> None:
        """Run parallel command"""
        parallel_cmd = self.get_parallel_cmd(
            cmd_list, parallel_cmds_file, parallel_log_file
        )
        sp.run(parallel_cmd, shell=True)
        Path(parallel_cmds_file).unlink()

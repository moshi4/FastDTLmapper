from pathlib import Path

import pytest
from fastdtlmapper.args import Args
from fastdtlmapper.cmd import Cmd


@pytest.fixture()
def args():
    """args fixture"""
    return Args(
        indir=Path("./testdir/"),
        tree_file=Path("./path/tree.nwk"),
        outdir=Path("./outdir/"),
        process_num=1,
        dup_cost=2,
        los_cost=1,
        trn_cost=1,
        inflation=1.5,
        timetree=False,
        rseed=0,
    )


def test_get_make_ultrametric_cmd(args: Args):
    """test get_make_ultrametric_cmd"""
    nwk_tree_infile, nwk_tree_outfile = "test_infile.nwk", "test_outfile.nwk"
    cmd = Cmd(args).get_make_ultrametric_cmd(nwk_tree_infile, nwk_tree_outfile)
    expected_cmd = f"make_ultrametric.py {nwk_tree_infile} {nwk_tree_outfile} "
    assert cmd == expected_cmd


def test_get_orthofinder_cmd(args: Args):
    """test get_orthofinder_cmd"""
    fasta_indir = "fasta_indir"
    cmd = Cmd(args).get_orthofinder_cmd(fasta_indir)
    expected_cmd = (
        f"orthofinder.py -og -f {fasta_indir} -t {args.process_num} "
        + f"-I {args.inflation}"
    )
    assert cmd == expected_cmd


def test_get_mafft_cmd(args: Args):
    """test get_mafft_cmd"""
    fasta_infile, aln_outfile = "infile.fa", "outfile_aln.fa"
    cmd = Cmd(args).get_mafft_cmd(fasta_infile, aln_outfile)
    expected_cmd = (
        f"mafft --auto --anysymbol --quiet {fasta_infile} > {aln_outfile} 2>&1"
    )
    assert cmd == expected_cmd


def test_get_trimal_cmd(args: Args):
    """test get_trimal_cmd"""
    aln_infile, aln_trim_outfile = "infile_aln.fa", "outfile_aln_trim.fa"
    cmd = Cmd(args).get_trimal_cmd(aln_infile, aln_trim_outfile)
    expected_cmd = f"trimal -in {aln_infile} -out {aln_trim_outfile} -automated1 2>&1"
    return cmd == expected_cmd


def test_get_iqtree_cmd_boot_option_on(args: Args):
    """test get_iqtree_cmd boot option on"""
    aln_infile, prefix = "infile_aln.fa", "prefix"
    cmd = Cmd(args).get_iqtree_cmd(aln_infile, prefix, True)
    expected_cmd = (
        f"iqtree -s {aln_infile} --prefix {prefix} -m TEST -mset JTT,WAG,LG "
        + f"--seed {args.rseed} --ufboot 1000 --boot-trees --wbtl --redo --quiet 2>&1"
    )
    assert cmd == expected_cmd


def test_get_iqtree_cmd_boot_option_off(args: Args):
    """test get_iqtree_cmd boot option off"""
    aln_infile, prefix = "infile_aln.fa", "prefix"
    cmd = Cmd(args).get_iqtree_cmd(aln_infile, prefix, False)
    expected_cmd = (
        f"iqtree -s {aln_infile} --prefix {prefix} -m TEST -mset JTT,WAG,LG "
        + f"--seed {args.rseed}  --redo --quiet 2>&1"
    )
    assert cmd == expected_cmd


def test_get_angst_cmd_timetree_option_on(args: Args):
    """Get AnGST run command"""
    species_tree_file, boot_tree_file, outdir = "tree.nwk", "boottree.nwk", "outdir"
    args.timetree = True
    cmd = Cmd(args).get_angst_cmd(species_tree_file, boot_tree_file, outdir)
    expected_cmd = (
        f"AnGST_wrapper.py -s {species_tree_file} -g {boot_tree_file} -o {outdir} "
        + f"--dup_cost {args.dup_cost} --los_cost {args.los_cost} "
        + f"--trn_cost {args.trn_cost} --timetree 2>&1"
    )
    assert cmd == expected_cmd


def test_get_parallel_cmd(args: Args, tmp_path: Path):
    """test get_parallle_cmd"""
    cmd_list = ["echo 1", "echo 2", "echo 3"]
    testdir = tmp_path / "test"
    testdir.mkdir()
    parallel_cmds_file = testdir / "parallel_cmds.txt"
    parallel_log_file = testdir / "parallel_log.txt"
    cmd = Cmd(args).get_parallel_cmd(cmd_list, parallel_cmds_file, parallel_log_file)
    expected_cmd = (
        f"parallel --no-notice --bar -a {parallel_cmds_file} "
        + f"-j {args.process_num} --results {parallel_log_file} "
        + "> /dev/null "
    )
    assert cmd == expected_cmd and parallel_cmds_file.is_file()


def test_run_parallel_cmd(args: Args, tmp_path: Path):
    """test run_parallel_cmd"""
    cmd_list = ["echo 1", "echo 2", "echo 3"]
    testdir = tmp_path / "test"
    testdir.mkdir()
    parallel_cmds_file = testdir / "parallel_cmds.txt"
    parallel_log_file = testdir / "parallel_log.csv"
    Cmd(args).run_parallel_cmd(cmd_list, parallel_cmds_file, parallel_log_file)
    assert not parallel_cmds_file.is_file() and parallel_log_file.is_file()

import os
from pathlib import Path

from fastdtlmapper.args import get_args


def test_get_args_default_ok():
    """OK Case: default parameters"""
    indir, tree_file, outdir = "indir/seq/", "indir/tree.nwk", "outdir/"
    argv = f"-i {indir} -t {tree_file} -o {outdir}"
    args = get_args(argv.split(" "))
    cpu_count = os.cpu_count()
    default_processor_num = cpu_count - 1 if cpu_count is not None else 1
    assert args.indir == Path(indir)
    assert args.tree_file == Path(tree_file)
    assert args.outdir == Path(outdir)
    assert args.process_num == default_processor_num
    assert args.dup_cost == 2
    assert args.los_cost == 1
    assert args.trn_cost == 3
    assert args.inflation == 3.0
    assert args.timetree is False
    assert args.rseed == 0


def test_get_args_user_specified_ok():
    """OK Case: all parameter user specified"""
    indir, tree_file, outdir = "indir/seq/", "indir/tree.nwk", "outdir/"
    process_num, dup_cost, los_cost, trn_cost = 4, 4, 2, 6
    inflation, rseed = 6.0, 100
    argv = (
        f"-i {indir} -t {tree_file} -o {outdir} -p {process_num} "
        + f"--dup_cost {dup_cost} --los_cost {los_cost} --trn_cost {trn_cost} "
        + f"--inflation {inflation} --timetree --rseed {rseed}"
    )
    args = get_args(argv.split(" "))
    assert args.indir == Path(indir)
    assert args.tree_file == Path(tree_file)
    assert args.outdir == Path(outdir)
    assert args.process_num == process_num
    assert args.dup_cost == dup_cost
    assert args.los_cost == los_cost
    assert args.trn_cost == trn_cost
    assert args.inflation == inflation
    assert args.timetree is True
    assert args.rseed == 100

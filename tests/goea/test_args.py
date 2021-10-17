from pathlib import Path

import pytest
from fastdtlmapper.goea.args import get_args


def test_get_args_default_ok():
    """OK Case: default parameters"""
    indir = "inputdir/"
    argv = f"-i {indir}"
    args = get_args(argv.split(" "))
    assert args.indir == Path(indir)
    assert args.plot_pvalue_thr == 0.05
    assert args.plot_max_num == 10
    assert args.plot_format == "png"
    assert args.plot_color == ""
    assert args.use_adjusted_pvalue is False


def test_get_args_user_specified_ok():
    """OK Case: all parameter user specified"""
    indir = "inputdir/"
    plot_pvalue_thr = 0.01
    plot_max_num = 20
    plot_format = "jpg"
    plot_color = "1affdb"

    argv = (
        f"-i {indir} --plot_pvalue_thr {plot_pvalue_thr} "
        + f"--plot_max_num {plot_max_num} --plot_format {plot_format} "
        + f"--plot_color {plot_color} --adjusted_pvalue"
    )
    args = get_args(argv.split(" "))
    assert args.indir == Path(indir)
    assert args.plot_pvalue_thr == plot_pvalue_thr
    assert args.plot_max_num == plot_max_num
    assert args.plot_format == plot_format
    assert args.plot_color == f"#{plot_color}"
    assert args.use_adjusted_pvalue is True


def test_get_args_invalid_plot_color_check():
    """NG case: Incorrect plot color specified"""
    indir = "inputdir/"
    # Check01. More than 6 characters
    # Check02. Contain no color code (z)
    # Check03. Less than 6 characters
    invalid_plot_color_list = ["de11s3a", "fei3zz", "de3a"]
    for invalid_plot_color in invalid_plot_color_list:
        argv = f"-i {indir} --plot_color {invalid_plot_color}"
        with pytest.raises(SystemExit) as e:
            get_args(argv.split(" "))
        assert e.type == SystemExit
        assert e.value.code == 2

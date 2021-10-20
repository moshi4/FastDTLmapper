from pathlib import Path

import pytest
from fastdtlmapper.scripts.plot_gain_loss_map import get_args, main


def test_plot_gain_loss_map_default_param(
    gain_loss_map_nwk_file: Path,
    tmp_path: Path,
):
    """test plot_gain_loss_map with default param"""
    outfile = tmp_path / "all_gain_loss_map.png"
    argv = f"-i {gain_loss_map_nwk_file} -o {outfile}"
    args = get_args(argv.split(" "))
    main(args)
    assert outfile.is_file()


def test_plot_gain_loss_map_user_param(
    gain_loss_map_nwk_file: Path,
    tmp_path: Path,
):
    """test plot_gain_loss_map with user param"""
    outfile = tmp_path / "all_gain_loss_map.svg"
    # Perfect plot test required confirming figure contents
    # It's impossible, so Color,Symbol,FontSize options are ignored
    argv = (
        f"-i {gain_loss_map_nwk_file} -o {outfile} --plot_margin 100 "
        + "--plot_width 1000 --title helloworld --ladderize"
    )
    args = get_args(argv.split(" "))
    main(args)
    assert outfile.is_file()


def test_plot_gain_loss_map_ng_extension(
    gain_loss_map_nwk_file: Path,
    tmp_path: Path,
):
    """test plot_gain_loss_map with outfile extension error"""
    outfile_invalid_ext = tmp_path / "all_gain_loss_map.invalid"
    argv = f"-i {gain_loss_map_nwk_file} -o {outfile_invalid_ext}"
    with pytest.raises(SystemExit) as e:
        get_args(argv.split(" "))
    assert e.type == SystemExit
    assert e.value.code == 2

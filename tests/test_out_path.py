from pathlib import Path

import pytest
from fastdtlmapper.out_path import OutPath


def test_makedirs(tmp_path: Path):
    """test OutPath makedirs"""
    outpath = OutPath(tmp_path)
    # Check makedirs results
    assert outpath.user_fasta_dir.is_dir()
    assert outpath.user_fasta_dir.is_dir()
    assert outpath.ortho_dir.is_dir()
    assert outpath.dtl_rec_dir.is_dir()
    assert outpath.aggregate_map_dir.is_dir()
    assert outpath.parallel_cmds_dir.is_dir()


def test_makedirs_goea_ok(tmp_path: Path):
    """test OutPath makedirs GOEA mode OK case"""
    # Create FastDTLmapper result directory
    outpath = OutPath(tmp_path)
    # No error occur, after FastDTLmapper run
    OutPath(outpath.rootdir, goea_mode=True)
    assert outpath.go_annotation_workdir.is_dir()
    assert outpath.go_enrichment_dir.is_dir()
    assert outpath.result_summary_plot_dir.is_dir()


def test_makedirs_goea_ng(tmp_path: Path):
    """test OutPath makedirs GOEA mode NG case"""
    with pytest.raises(ValueError):
        # Raised error if user specifiy no FastDTLmapper result dir
        OutPath(tmp_path, goea_mode=True)

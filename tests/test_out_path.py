from pathlib import Path

from fastdtlmapper.out_path import OutPath


def test_makedirs(tmp_path: Path):
    """test OutPath makedirs"""
    rootdir = tmp_path / "rootdir"
    rootdir.mkdir()
    outpath = OutPath(rootdir)
    assert (
        outpath.user_fasta_dir.is_dir()
        and outpath.user_fasta_dir.is_dir()
        and outpath.ortho_dir.is_dir()
        and outpath.dtl_rec_dir.is_dir()
        and outpath.aggregate_map_dir.is_dir()
        and outpath.parallel_cmds_dir.is_dir()
    )

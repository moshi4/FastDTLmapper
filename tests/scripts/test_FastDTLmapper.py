from pathlib import Path

from fastdtlmapper.args import get_args
from fastdtlmapper.scripts import FastDTLmapper


def test_integration_fastdtlmapper_main(
    integration_fasta_indir: Path,
    integration_species_tree_file: Path,
    tmp_path: Path,
):
    """FastDTLmapper integration run test"""
    argv = (
        f"-i {integration_fasta_indir} "
        + f"-t {integration_species_tree_file} "
        + f"-o {tmp_path}"
    )
    args = get_args(argv.split(" "))
    FastDTLmapper.main(args)

from pathlib import Path

import pytest
from fastdtlmapper.input_check import InputCheck


def test_input_check_run(fasta_or_genbank_dir: Path, species_tree_file: Path):
    """test input_check run ok"""
    InputCheck(fasta_or_genbank_dir, species_tree_file).run()


def test_get_fasta_or_genbank_file_list(
    fasta_or_genbank_dir: Path, species_tree_file: Path
):
    """test get_fasta_or_genbank_file_list count"""
    input_check = InputCheck(fasta_or_genbank_dir, species_tree_file)
    assert len(input_check._get_fasta_or_genbank_file_list()) == 7


def test_is_valid_filename_ng():
    """test is_valid_filename NG case"""
    input_check = InputCheck(Path("dummy_dir"), Path("dummy.nwk"))
    invalid_filename_list = [
        Path("species_test1.fa"),
        Path("species-test2.gbk"),
        Path("species|test3.faa"),
    ]
    with pytest.raises(SystemExit) as e:
        input_check._is_valid_filename(invalid_filename_list)
    assert e.type == SystemExit
    assert e.value.code == 1


def test_check_rooted_tree_ng(
    fasta_or_genbank_dir: Path, unrooted_species_tree_file: Path
):
    """test get_fasta_or_genbank_file_list count"""
    input_check = InputCheck(fasta_or_genbank_dir, unrooted_species_tree_file)
    with pytest.raises(SystemExit) as e:
        input_check._check_rooted_tree()
    assert e.type == SystemExit
    assert e.value.code == 1


def test_check_tree_seqfile_consistency_ng_diff_num(
    fasta_or_genbank_dir: Path, species_tree_file: Path
):
    """test check_tree_seqfile_consistency NG case (diff number)"""
    input_check = InputCheck(fasta_or_genbank_dir, species_tree_file)
    # Number of seqfile and tree species number different case
    with pytest.raises(SystemExit) as e:
        input_check._check_tree_seqfile_consistency([])
    assert e.type == SystemExit
    assert e.value.code == 1


def test_check_tree_seqfile_consistency_ng_name_unmatch(
    fasta_or_genbank_dir: Path, species_tree_file: Path
):
    """test check_tree_seqfile_consistency NG case (diff number)"""
    input_check = InputCheck(fasta_or_genbank_dir, species_tree_file)
    # Seqfile names & species names unmatch case
    with pytest.raises(SystemExit) as e:
        input_check._check_tree_seqfile_consistency([Path("test")] * 7)
    assert e.type == SystemExit
    assert e.value.code == 1

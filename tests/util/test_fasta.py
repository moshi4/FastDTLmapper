from pathlib import Path

import pytest
from fastdtlmapper.util import UtilFasta
from pytest import TempPathFactory


@pytest.fixture()
def fasta_file(tmp_path_factory: TempPathFactory) -> Path:
    """Fasta file"""
    fasta_text = "\n".join(
        [
            ">sp1_GENE000001",
            "X" * 50,
            ">sp2_GENE000002",
            "X" * 100,
            ">sp3_GENE000003",
            "X" * 100,
            ">sp3_GENE000003",
            "Y" * 100,
        ]
    )
    fasta_file_path = tmp_path_factory.mktemp("fasta") / "seq.fa"
    fasta_file_path.write_text(fasta_text)
    return fasta_file_path


@pytest.fixture()
def blank_file(tmp_path_factory: TempPathFactory) -> Path:
    """Blank file"""
    blank_file_path = tmp_path_factory.mktemp("fasta") / "blank.txt"
    blank_file_path.write_text("")
    return blank_file_path


def test_seq_num(fasta_file: Path) -> None:
    """test seq_num"""
    assert UtilFasta(fasta_file).seq_num == 4


def test_uniq_seq_num(fasta_file: Path) -> None:
    """test uniq_seq_num"""
    assert UtilFasta(fasta_file).uniq_seq_num == 3


def test_is_valid_format_ok(fasta_file: Path) -> None:
    """test is_valid_format ok"""
    assert UtilFasta(fasta_file).is_valid_format


def test_is_valid_format_ng(blank_file: Path) -> None:
    """test is_valid_format ng"""
    assert not UtilFasta(blank_file).is_valid_format


def test_species_name2seq_num(fasta_file: Path) -> None:
    """test species_name2seq_num"""
    assert UtilFasta(fasta_file).species_name2seq_num == {
        "sp1": 1,
        "sp2": 1,
        "sp3": 2,
    }


def test_uniq_species_name_list(fasta_file: Path) -> None:
    """test uniq_species_name_list"""
    assert UtilFasta(fasta_file).uniq_species_name_list == ["sp1", "sp2", "sp3"]

from pathlib import Path

from fastdtlmapper.util import UtilFasta


def test_seq_num(fasta_file: Path):
    """test seq_num"""
    assert UtilFasta(fasta_file).seq_num == 4


def test_uniq_seq_num(fasta_file: Path):
    """test uniq_seq_num"""
    assert UtilFasta(fasta_file).uniq_seq_num == 3


def test_is_valid_format_ok(fasta_file: Path):
    """test is_valid_format ok"""
    assert UtilFasta(fasta_file).is_valid_format


def test_is_valid_format_ng(blank_file: Path):
    """test is_valid_format ng"""
    assert not UtilFasta(blank_file).is_valid_format


def test_species_name2seq_num(fasta_file: Path):
    """test species_name2seq_num"""
    assert UtilFasta(fasta_file).species_name2seq_num == {
        "sp1": 1,
        "sp2": 1,
        "sp3": 2,
    }


def test_uniq_species_name_list(fasta_file: Path):
    """test uniq_species_name_list"""
    assert UtilFasta(fasta_file).uniq_species_name_list == ["sp1", "sp2", "sp3"]


def test_add_serial_id(fasta_file: Path, tmp_path: Path):
    """test add_serial_id"""
    fasta_serial_id_outfile = tmp_path / "serial_id.fa"
    UtilFasta(fasta_file).add_serial_id(fasta_serial_id_outfile, "prefix")
    assert fasta_serial_id_outfile.is_file()

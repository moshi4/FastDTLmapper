from pathlib import Path

from Bio import SeqIO
from fastdtlmapper.util import UtilGenbank


def test_is_valid_format_ok(genbank_file: Path):
    """test is_valid_format OK case"""
    assert UtilGenbank(genbank_file).is_valid_format


def test_is_valid_format_ng(blank_file: Path):
    """test is_valid_format NG case"""
    assert not UtilGenbank(blank_file).is_valid_format


def test_convert_cds_fasta_protein(genbank_file: Path, tmp_path: Path):
    """test convert cds_fasta (seqtype=protein)"""
    cds_protein_fasta_outfile = tmp_path / "cds.fa"
    UtilGenbank(genbank_file).convert_cds_fasta(
        cds_protein_fasta_outfile, "protein", "prefix"
    )
    assert len(list(SeqIO.parse(cds_protein_fasta_outfile, "fasta"))) != 0


def test_convert_cds_fasta_nucleotide(genbank_file: Path, tmp_path: Path):
    """test convert cds_fasta (seqtype=nucleotide)"""
    cds_nucleotide_fasta_outfile = tmp_path / "cds.fa"
    UtilGenbank(genbank_file).convert_cds_fasta(
        cds_nucleotide_fasta_outfile, "nucleotide", "prefix"
    )
    assert len(list(SeqIO.parse(cds_nucleotide_fasta_outfile, "fasta"))) != 0

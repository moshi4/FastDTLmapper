from pathlib import Path

import pytest
from fastdtlmapper.angst.model.transfer import Transfer
from fastdtlmapper.reconcilation import Reconciliation as Rec


@pytest.fixture()
def species_tree_file(tmp_path: Path) -> Path:
    """species tree file fixture"""
    file = tmp_path / "species_tree.nwk"
    tree_text = (
        "((outsp1,outsp2)N002,"
        + "((insp1,insp2)N004,(insp3,(insp4,insp5)N006)N005)N003)N001;"
    )
    file.write_text(tree_text)
    return file


@pytest.fixture()
def fasta_2seq_file_insp1_only(tmp_path: Path) -> Path:
    """2seq 1species fasta file fixture"""
    file = tmp_path / "seq.fa"
    fasta_text = f">insp1_GENE000001\n{'X'*100}\n" + f">insp1_GENE000002\n{'X'*100}"
    file.write_text(fasta_text)
    return file


def test_two_species_brn_dup(species_tree_file: Path, fasta_2seq_file_insp1_only: Path):
    """test one species brn & dup case"""
    node_event_list = Rec.two_species(species_tree_file, fasta_2seq_file_insp1_only)
    for node_event in node_event_list:
        if node_event.node_id == "insp1":
            assert (
                node_event.brn_num == 1
                and node_event.dup_num == 1
                and node_event.los_num == 0
                and node_event.trn_num == 0
                and node_event.trn_detail_list == []
                and node_event.gene_num == 2
            )
        else:
            assert (
                node_event.brn_num == 0
                and node_event.dup_num == 0
                and node_event.los_num == 0
                and node_event.trn_num == 0
                and node_event.trn_detail_list == []
                and node_event.gene_num == 0
            )


@pytest.fixture()
def fasta_2seq_file_insp1_insp3(tmp_path: Path) -> Path:
    """2seq 1species fasta file fixture"""
    file = tmp_path / "seq.fa"
    fasta_text = f">insp1_GENE000001\n{'X'*100}\n" + f">insp3_GENE000002\n{'X'*100}"
    file.write_text(fasta_text)
    return file


def test_two_species_brn_los(
    species_tree_file: Path, fasta_2seq_file_insp1_insp3: Path
):
    """test two species brn & los case"""
    node_event_list = Rec.two_species(species_tree_file, fasta_2seq_file_insp1_insp3)
    for node_event in node_event_list:
        if node_event.node_id == "N003":
            assert (
                node_event.brn_num == 1
                and node_event.dup_num == 0
                and node_event.los_num == 0
                and node_event.trn_num == 0
                and node_event.trn_detail_list == []
                and node_event.gene_num == 1
            )

        elif node_event.node_id in ("insp1", "insp3", "N004", "N005"):
            assert (
                node_event.brn_num == 0
                and node_event.dup_num == 0
                and node_event.los_num == 0
                and node_event.trn_num == 0
                and node_event.trn_detail_list == []
                and node_event.gene_num == 1
            )
        elif node_event.node_id in ("insp2", "N006"):
            assert (
                node_event.brn_num == 0
                and node_event.dup_num == 0
                and node_event.los_num == 1
                and node_event.trn_num == 0
                and node_event.trn_detail_list == []
                and node_event.gene_num == 0
            )
        else:
            assert (
                node_event.brn_num == 0
                and node_event.dup_num == 0
                and node_event.los_num == 0
                and node_event.trn_num == 0
                and node_event.trn_detail_list == []
                and node_event.gene_num == 0
            )


@pytest.fixture()
def fasta_2seq_file_insp5_outsp1(tmp_path: Path) -> Path:
    """2seq 1species fasta file fixture"""
    file = tmp_path / "seq.fa"
    fasta_text = f">insp5_GENE000001\n{'X'*100}\n" + f">outsp1_GENE000002\n{'X'*100}"
    file.write_text(fasta_text)
    return file


def test_two_species_brn_trn(
    species_tree_file: Path, fasta_2seq_file_insp5_outsp1: Path
):
    """test two species brn & los case"""
    node_event_list = Rec.two_species(species_tree_file, fasta_2seq_file_insp5_outsp1)
    for node_event in node_event_list:
        if node_event.node_id == "insp5":
            assert (
                node_event.brn_num == 1
                and node_event.dup_num == 0
                and node_event.los_num == 0
                and node_event.trn_num == 0
                and node_event.trn_detail_list == []
                and node_event.gene_num == 1
            )
        elif node_event.node_id == "outsp1":
            assert (
                node_event.brn_num == 0
                and node_event.dup_num == 0
                and node_event.los_num == 0
                and node_event.trn_num == 1
                and node_event.trn_detail_list == [Transfer("insp5", "outsp1", False)]
                and node_event.gene_num == 1
            )
        else:
            assert (
                node_event.brn_num == 0
                and node_event.dup_num == 0
                and node_event.los_num == 0
                and node_event.trn_num == 0
                and node_event.trn_detail_list == []
                and node_event.gene_num == 0
            )


@pytest.fixture()
def fasta_1seq_file(tmp_path: Path) -> Path:
    """1seq fasta file fixture"""
    file = tmp_path / "seq.fa"
    fasta_text = f">insp1_GENE000001\n{'X'*100}"
    file.write_text(fasta_text)
    return file


def test_one_species(species_tree_file: Path, fasta_1seq_file: Path):
    """test one_species"""
    node_event_list = Rec.one_species(species_tree_file, fasta_1seq_file)
    for node_event in node_event_list:
        if node_event.node_id == "insp1":
            assert (
                node_event.brn_num == 1
                and node_event.dup_num == 0
                and node_event.los_num == 0
                and node_event.trn_num == 0
                and node_event.trn_detail_list == []
                and node_event.gene_num == 1
            )
        else:
            assert (
                node_event.brn_num == 0
                and node_event.dup_num == 0
                and node_event.los_num == 0
                and node_event.trn_num == 0
                and node_event.trn_detail_list == []
                and node_event.gene_num == 0
            )

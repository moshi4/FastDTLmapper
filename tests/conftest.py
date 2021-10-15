from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def data_dir() -> Path:
    """data directory fixture"""
    return Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def genbank_file(data_dir: Path) -> Path:
    """genbank file fixture"""
    return data_dir / "test.gbk"


@pytest.fixture(scope="session")
def fasta_file(data_dir: Path) -> Path:
    """fasta file fixture"""
    return data_dir / "test.fa"


@pytest.fixture(scope="session")
def blank_file(data_dir: Path) -> Path:
    """blank file fixture"""
    return data_dir / "blank.txt"


@pytest.fixture(scope="session")
def fasta_or_genbank_dir(data_dir: Path) -> Path:
    """fasta or genbank directory fixture"""
    return data_dir / "fasta_or_genbank"


@pytest.fixture(scope="session")
def species_tree_file(data_dir: Path) -> Path:
    """species tree file fixture"""
    return data_dir / "tree" / "species_tree.nwk"


@pytest.fixture(scope="session")
def unrooted_species_tree_file(data_dir: Path) -> Path:
    """unrooted species tree file fixture"""
    return data_dir / "tree" / "unrooted_species_tree.nwk"


@pytest.fixture(scope="session")
def species_tree_nodeid_file(data_dir: Path) -> Path:
    """species tree with nodeid file fixture"""
    return data_dir / "tree" / "species_tree_nodeid.nwk"


@pytest.fixture(scope="session")
def angst_result_dir(data_dir: Path) -> Path:
    """AnGST output directory fixture"""
    return data_dir / "angst"

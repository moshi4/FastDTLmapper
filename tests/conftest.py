from pathlib import Path
from typing import List

import pytest


@pytest.fixture(scope="session")
def data_dir() -> Path:
    """data directory fixture"""
    return Path(__file__).parent / "data"


###########################################################
# FastDTLmapper testdata fixture
###########################################################


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


@pytest.fixture(scope="session")
def integration_fasta_indir(data_dir: Path) -> Path:
    """FastDTLmapper integration fasta dir"""
    return data_dir / "integration_test" / "fasta"


@pytest.fixture(scope="session")
def integration_species_tree_file(data_dir: Path) -> Path:
    """FastDTLmapper integration species tree file"""
    return data_dir / "integration_test" / "species_tree.nwk"


###############################################################################
# FastDTLgoea testdata fixture
###############################################################################


@pytest.fixture(scope="session")
def orthogroups_file(data_dir: Path) -> Path:
    """OrthoFinder ortholog groups file fixture"""
    return data_dir / "goea" / "Orthogroups.txt"


@pytest.fixture(scope="session")
def go_annotation_file_list(data_dir: Path) -> List[Path]:
    """Interproscan go annotation file list fixture"""
    go_annotation_dir = data_dir / "goea" / "go_annotation"
    return list(go_annotation_dir.glob("*.tsv"))


@pytest.fixture(scope="session")
def go_basic_obo_file(data_dir: Path) -> Path:
    """go-basic.obo file fixture"""
    return data_dir / "goea" / "go-basic.obo"


@pytest.fixture(scope="session")
def goea_target_gene_list_file(data_dir: Path) -> Path:
    """GOEA target gene list file fixture"""
    return data_dir / "goea" / "goea_target_gene_list.txt"


@pytest.fixture(scope="session")
def goea_all_gene_list_file(data_dir: Path) -> Path:
    """GOEA all gene list file fixture"""
    return data_dir / "goea" / "goea_all_gene_list.txt"


@pytest.fixture(scope="session")
def goea_association_file(data_dir: Path) -> Path:
    """GOEA association file fixture"""
    return data_dir / "goea" / "goea_association.txt"


@pytest.fixture(scope="session")
def goea_result_file(data_dir: Path) -> Path:
    """GOEA goatools result file fixture"""
    return data_dir / "goea" / "goea_result.tsv"


###########################################################
# plot_gain_loss_map testdata fixture
###########################################################


@pytest.fixture(scope="session")
def dtl_map_nwk_file(data_dir: Path) -> Path:
    """dtl map newick file fixture"""
    return data_dir / "all_dtl_map.nwk"

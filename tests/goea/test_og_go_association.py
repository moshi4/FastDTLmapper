from pathlib import Path
from typing import Dict, List, Set

import pytest
from fastdtlmapper.goea.og_go_association import OgGoAssociation


def test_get_og_id2gene_id_association(tmp_path: Path):
    """test get_og_id2gene_id_association"""
    # Dummy test input orthogroups file
    og_gene_list_file = tmp_path / "Orthogroups.txt"
    og_gene_list_file.write_text(
        "OG0000000: GENE000001 GENE000003 GENE000005\n"
        "OG0000001: GENE000007 GENE000008 GENE000009 GENE000010\n"
    )
    # Check orthogroups file parse result
    oga = OgGoAssociation(og_gene_list_file, [])
    gene_id2go_id_association = oga._get_og_id2gene_id_association()
    assert gene_id2go_id_association["OG0000000"] == [
        "GENE000001",
        "GENE000003",
        "GENE000005",
    ]
    assert gene_id2go_id_association["OG0000001"] == [
        "GENE000007",
        "GENE000008",
        "GENE000009",
        "GENE000010",
    ]


def test_get_gene_id2go_id_association(tmp_path: Path):
    """test get_gene_id2go_id_association"""
    # Dummy test input gene annotation file
    gene_annotation_file = tmp_path / "gene_annotation.tsv"
    gene_annotation_file.write_text(
        ("GENE000001\t" + "-\t" * 12 + "-\n")
        + ("GENE000001\t" + "-\t" * 12 + "GO:0000001|GO:0000002\n")
        + ("GENE000001\t" + "-\t" * 12 + "GO:0000003\n")
        + ("GENE000001\t" + "-\t" * 12 + "GO:0000002|GO:0000003\n")
        + ("GENE000002\t" + "-\t" * 11 + "\n")
        + ("GENE000002\t" + "-\t" * 12 + "GO:0000001|GO:0000004\n")
    )
    # Check gene annotation file parse result
    oga = OgGoAssociation(Path("dummy"), [gene_annotation_file])
    gene_id2go_id_association = oga._get_gene_id2go_id_association()
    assert sorted(gene_id2go_id_association["GENE000001"]) == sorted(
        set(["GO:0000001", "GO:0000002", "GO:0000003"])
    )
    assert sorted(gene_id2go_id_association["GENE000002"]) == sorted(
        set(["GO:0000001", "GO:0000004"])
    )


@pytest.fixture()
def og_id2gene_id_list() -> Dict[str, List[str]]:
    """og_id2gene_id_list testdata fixture"""
    return {
        "OG0000000": ["GENE000001", "GENE000002", "GENE000003"],
        "OG0000001": ["GENE000004", "GENE000005"],
    }


@pytest.fixture()
def gene_id2go_id_set() -> Dict[str, Set[str]]:
    """gene_id2go_id_set testdata fixture"""
    return {
        "GENE000001": set(["GO:0000001", "GO:0000002", "GO:0000004"]),
        "GENE000002": set(["GO:0000001", "GO:0000002", "GO:0000003", "GO:0000004"]),
        "GENE000003": set(["GO:0000001", "GO:0000002", "GO:0000005"]),
        "GENE000004": set(["GO:0000001", "GO:0000002"]),
        "GENE000005": set(["GO:0000001", "GO:0000003", "GO:0000005"]),
    }


def test_og_id2go_id_association_ratio_half(
    og_id2gene_id_list: Dict[str, List[str]],
    gene_id2go_id_set: Dict[str, Set[str]],
    go_define_ratio_thr=0.5,
):
    """test get_og_id2go_id_association (ratio=0.5 [optimal threshold])"""
    oga = OgGoAssociation(Path("dummy"), [], go_define_ratio_thr)
    og_id2go_id_association = oga._get_og_id2go_id_association(
        og_id2gene_id_list, gene_id2go_id_set
    )
    assert sorted(og_id2go_id_association["OG0000000"]) == sorted(
        ["GO:0000001", "GO:0000002", "GO:0000004"]
    )
    assert sorted(og_id2go_id_association["OG0000001"]) == sorted(
        ["GO:0000001", "GO:0000002", "GO:0000003", "GO:0000005"]
    )


def test_og_id2go_id_association_ratio_one(
    og_id2gene_id_list: Dict[str, List[str]],
    gene_id2go_id_set: Dict[str, Set[str]],
    go_define_ratio_thr=1.0,
):
    """test get_og_id2go_id_association (ratio=1.0 [strict threshold])"""
    oga = OgGoAssociation(Path("dummy"), [], go_define_ratio_thr)
    og_id2go_id_association = oga._get_og_id2go_id_association(
        og_id2gene_id_list, gene_id2go_id_set
    )
    assert sorted(og_id2go_id_association["OG0000000"]) == sorted(
        ["GO:0000001", "GO:0000002"]
    )
    assert sorted(og_id2go_id_association["OG0000001"]) == sorted(
        ["GO:0000001"],
    )


def test_og_id2go_id_association_ratio_zero(
    og_id2gene_id_list: Dict[str, List[str]],
    gene_id2go_id_set: Dict[str, Set[str]],
    go_define_ratio_thr=0.0,
):
    """test get_og_id2go_id_association (ratio=0.0 [loose threshold])"""
    oga = OgGoAssociation(Path("dummy"), [], go_define_ratio_thr)
    og_id2go_id_association = oga._get_og_id2go_id_association(
        og_id2gene_id_list, gene_id2go_id_set
    )
    assert sorted(og_id2go_id_association["OG0000000"]) == sorted(
        ["GO:0000001", "GO:0000002", "GO:0000003", "GO:0000004", "GO:0000005"]
    )
    assert sorted(og_id2go_id_association["OG0000001"]) == sorted(
        ["GO:0000001", "GO:0000002", "GO:0000003", "GO:0000005"]
    )


def test_write_og2go_association(
    orthogroups_file: Path,
    go_annotation_file_list: List[Path],
    tmp_path: Path,
):
    """test write_og2go_association (Check run successfully)"""
    oga = OgGoAssociation(orthogroups_file, go_annotation_file_list)
    og2go_association_file = tmp_path / "og2go_association.txt"
    oga.write_og2go_association(og2go_association_file)
    assert og2go_association_file.is_file()

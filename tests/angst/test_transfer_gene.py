from pathlib import Path

from fastdtlmapper.angst import AngstTransferGene


def test_transfer_gene(species_tree_nodeid_file: Path, angst_result_dir: Path):
    """test angst_transfer_gene"""
    trn_gene = AngstTransferGene(species_tree_nodeid_file, angst_result_dir)
    assert len(trn_gene.trn_fromto2gene_id_list.keys()) == 1
    assert trn_gene.trn_fromto2gene_id_list["N002 -> insp3"] == [
        "insp3_GENE000634",
        "insp3_GENE000063",
    ]

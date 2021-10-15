from pathlib import Path

from fastdtlmapper.util import UtilTree


def test_is_valid_newick_format_ok(species_tree_file: Path):
    """test is_valid_newick_format OK case"""
    assert UtilTree(species_tree_file).is_valid_newick_format()


def test_is_valid_newick_format_ng(blank_file: Path):
    """test is_valid_newick_format NG case"""
    assert not UtilTree(blank_file).is_valid_newick_format()


def test_is_rooted_tree_ok(species_tree_file: Path):
    """test is_rooted_tree OK case"""
    assert UtilTree(species_tree_file).is_rooted_tree()


def test_is_rooted_tree_ng(unrooted_species_tree_file: Path):
    """test is_rooted_tree NGJ case"""
    assert not UtilTree(unrooted_species_tree_file).is_rooted_tree()


def test_add_internal_node_id(species_tree_file: Path, tmp_path: Path):
    """test add_internal_node_id"""
    species_tree_nodeid_file = tmp_path / "species_tree_nodeid.nwk"
    UtilTree.add_internal_node_id(species_tree_file, species_tree_nodeid_file)
    # Check internal node id exists
    outtree_text = species_tree_nodeid_file.read_text()
    for node_id in ("N001", "N002", "N003", "N004", "N005", "N006"):
        assert outtree_text.count(node_id) == 1


def test_make_3gene_tree(tmp_path: Path):
    """test make_3gene_tree"""
    fasta_3gene_file = tmp_path / "3gene.fa"
    fasta_3gene_file.write_text(">insp1\nXXXXX\n>insp2\nXXXXX\n>outsp1\nXXXXX\n")
    gene_tree_outfile = tmp_path / "3gene_tree.nwk"
    UtilTree.make_3genes_tree(fasta_3gene_file, gene_tree_outfile)
    gene_tree_text = gene_tree_outfile.read_text()
    assert gene_tree_text == "(insp1:0.01,insp2:0.01,outsp1:0.01);"

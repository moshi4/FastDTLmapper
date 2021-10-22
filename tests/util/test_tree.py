from pathlib import Path

from ete3 import Tree
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


def test_multifurcate_zero_length_nodes(tmp_path: Path):
    """test multifurcate_zero_length_nodes"""
    root_tree_file = tmp_path / "root_tree.nwk"
    root_tree_text = "(((a:0,b:0):0,c:0):1,(d:1,e:1):1);"
    root_tree_file.write_text(root_tree_text)
    root_tree = Tree(root_tree_text)

    multifurcate_outfile = tmp_path / "multifuracate_tree.nwk"
    UtilTree.multifurcate_zero_length_nodes(root_tree_file, multifurcate_outfile)
    multifurcate_tree = Tree(multifurcate_outfile.read_text())

    assert len(root_tree.get_common_ancestor("a", "b", "c").children) == 2
    assert len(multifurcate_tree.get_common_ancestor("a", "b", "c").children) == 3


def test_unroot_tree(tmp_path: Path):
    """test unroot_tree"""
    root_tree_file = tmp_path / "root_tree.nwk"
    root_tree_text = "(((a,b),c),(d,e));"
    root_tree_file.write_text(root_tree_text)

    unroot_tree_file = tmp_path / "unroot_tree.nwk"
    UtilTree.unroot_tree(root_tree_file, unroot_tree_file)
    unroot_tree = Tree(unroot_tree_file.read_text())

    assert len(Tree(root_tree_text).get_tree_root().children) == 2
    assert len(unroot_tree.children) == 3

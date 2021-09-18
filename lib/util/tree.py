from Bio import Phylo
from Bio.Phylo.BaseTree import Tree


def is_ultrametric_tree(nwk_tree_file: str) -> bool:
    """Check ultrametric tree or not

    Args:
        nwk_tree_file (str): Newick format tree file path

    Returns:
        bool: Ultrametric check result
    """
    tree: Tree = Phylo.read(nwk_tree_file, "newick")
    # Get all path (root -> leaf) branch length
    total_branch_length_list = []
    for leaf_node in tree.get_terminals():
        total_branch_length = 0
        for node in tree.get_path(leaf_node.name):
            if node.branch_length is None:
                return False
            total_branch_length += node.branch_length
        total_branch_length_list.append(total_branch_length)

    return len(set(total_branch_length_list)) == 1


def convert_bootstrap_float_to_int(nwk_tree_infile: str, nwk_tree_outfile: str) -> None:
    """Convert tree bootstrap value float to int (e.g. 0.978 -> 98)

    Args:
        nwk_tree_infile (str): Input newick format tree file path
        nwk_tree_outfile (str): Output newick format tree file path
    """
    tree: Tree = Phylo.read(nwk_tree_infile, "newick")
    for node in tree.get_nonterminals():
        if node.confidence is None or node.confidence > 1:
            continue
        node.confidence = int(float(node.confidence) * 100)
    Phylo.write(tree, nwk_tree_outfile, "newick", format_confidence="%d")

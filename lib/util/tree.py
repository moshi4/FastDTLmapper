from dataclasses import dataclass

from Bio import Phylo
from Bio.Phylo.BaseTree import Tree


@dataclass
class UtilTree:
    """Tree Utility Class"""

    @staticmethod
    def convert_bootstrap_float_to_int(
        nwk_tree_infile: str, nwk_tree_outfile: str
    ) -> None:
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

    @staticmethod
    def add_internal_node_id(nwk_tree_infile: str, nwk_tree_outfile: str) -> None:
        """Add internal node to serial id (e.g. N001,N002,...N0XX)

        Args:
            nwk_tree_infile (str): Input newick format tree file path
            nwk_tree_outfile (str): Output newick format tree file path
        """
        tree: Tree = Phylo.read(nwk_tree_infile, "newick")
        for cnt, node in enumerate(tree.get_nonterminals(), 1):
            node.name = f"N{cnt:03d}"
            node.confidence = None
        Phylo.write(tree, nwk_tree_outfile, "newick", format_branch_length="%1.6f")

        # Remove root unexpected branch length value
        with open(nwk_tree_outfile) as f:
            replace_tree_info = f.read().replace(":0.000000;", ";").replace("'", "")
        with open(nwk_tree_outfile, "w") as f:
            f.write(replace_tree_info)

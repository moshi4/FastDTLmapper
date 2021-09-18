from copy import deepcopy
from dataclasses import dataclass
from typing import Dict, Tuple

from Bio import Phylo
from Bio.Phylo.BaseTree import Tree


@dataclass
class MowgliTree:
    """Mowgli Tree Class"""

    species_tree_file: str

    def __post_init__(self):
        tree: Tree = Phylo.read(self.species_tree_file, "newick")
        pruned_tree = self._prune_outgroup_node(tree)
        self._numbered_tree, self.node_num2name = self._get_numbered_tree(pruned_tree)
        self._named_tree = self._convert_tree_node_name(
            self._numbered_tree, self.node_num2name
        )

    @property
    def numbered_tree(self) -> Tree:
        """Copied numbered tree"""
        return deepcopy(self._numbered_tree)

    @property
    def named_tree(self) -> Tree:
        """Copied named tree"""
        return deepcopy(self._named_tree)

    def _prune_outgroup_node(self, tree: Tree) -> Tree:
        """Prune Mowgli auto added outgroup

        Args:
            tree (Tree): Mowgli species tree

        Returns:
            Tree: Outgroup pruned tree
        """
        tree = deepcopy(tree)
        outgroup_clade = list(tree.find_clades({"name": "OUTGROUP_.*"}))[0]
        tree.prune(outgroup_clade)
        return tree

    def _get_numbered_tree(self, tree: Tree) -> Tuple[Tree, Dict[str, str]]:
        """Get Mowgli numbered node name tree (convert no numbered node name)

        Args:
            tree (Tree): Mowgli species tree

        Returns:
            Tuple[Tree, Dict[str, str]]: Numbered tree & convert dictionary
        """
        tree = deepcopy(tree)
        node_num2name = {}
        for node in tree.get_terminals():
            split_node_name = node.name.split("_")
            node_name = "_".join(split_node_name[0:-1])
            node_num = str(split_node_name[-1])
            node_num2name[node_num] = node_name
            node.name = node_num
        for cnt, node in enumerate(tree.get_nonterminals(), 1):
            node_num = str(node.confidence)
            node.confidence = None
            node_name = f"N{cnt:03d}"
            node_num2name[node_num] = node_name
            node.name = node_num
        return tree, node_num2name

    def _convert_tree_node_name(
        self, tree: Tree, node_num2name: Dict[str, str]
    ) -> Tree:
        """Convert tree numbered node to named node (e.g. '8' -> 'N003')

        Args:
            tree (Tree): Numbered name tree
            node_num2name (Dict[str, str]): number to name convert dictionary

        Returns:
            Tree: Node name converted tree
        """
        tree = deepcopy(tree)
        for node in tree.find_clades():
            node.name = node_num2name[node.name]
        return tree

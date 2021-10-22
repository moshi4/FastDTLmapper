from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Union

from Bio import Phylo, SeqIO
from Bio.Phylo.BaseTree import Tree
from ete3 import Tree as EteTree
from fastdtlmapper.angst.model import NodeEvent


@dataclass
class UtilTree:
    """Tree Utility Class"""

    tree_file: Union[str, Path]

    def is_valid_newick_format(self) -> bool:
        """Check newick format or not"""
        try:
            tree: Tree = Phylo.read(self.tree_file, "newick")
            tree_node_num = len(list(tree.find_clades()))
        except ValueError:
            return False
        return False if tree_node_num == 1 else True

    def is_rooted_tree(self) -> bool:
        """Check rooted tree or not"""
        if not self.is_valid_newick_format():
            return False
        tree: Tree = Phylo.read(self.tree_file, "newick")
        return True if tree.root.is_bifurcating() else False

    @staticmethod
    def add_internal_node_id(
        nwk_tree_infile: Union[str, Path],
        nwk_tree_outfile: Union[str, Path],
    ) -> None:
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

    @staticmethod
    def make_3genes_tree(
        fasta_file: Union[str, Path],
        gene_tree_file: Union[str, Path],
    ) -> None:
        """Make 3 genes unrooted newick tree file

        Args:
            fasta_file (str): Input fasta file path
            gene_tree_file (str): Output newick tree file path
        """
        gene_name_list = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            gene_name_list.append(record.id)

        gene_tree_info = "("
        for gene_name in gene_name_list:
            gene_tree_info += f"{gene_name}:0.01,"
        gene_tree_info = gene_tree_info.rstrip(",")
        gene_tree_info += ");"
        with open(gene_tree_file, "w") as f:
            f.write(gene_tree_info)

    @staticmethod
    def map_node_event(
        species_tree_infile: Union[str, Path],
        node_id2node_event: Dict[str, NodeEvent],
        nwk_tree_outfile: Union[str, Path],
        map_type: str,
    ):
        """Mapping node event to newick tree

        Args:
            species_tree_infile (str): Input species tree file path
            node_id2node_event (Dict[str, NodeEvent]): node id & node event dict
            nwk_tree_outfile (str): Output newick file path
            map_type (str): mapping type ("gain-loss" or "dtl")
        """
        tree: Tree = Phylo.read(species_tree_infile, "newick")
        for node in tree.find_clades():
            node_event = node_id2node_event[node.name]
            if map_type == "gain-loss":
                node.name = node.name + " | " + node_event.as_gain_loss_text
            elif map_type == "dtl":
                node.name = node.name + " | " + node_event.as_dtl_text
            else:
                raise ValueError("type must be 'gain-loss' or 'dtl'")

        Phylo.write(tree, nwk_tree_outfile, "newick", format_branch_length="%1.6f")

        # Remove root unexpected branch length value
        with open(nwk_tree_outfile) as f:
            replace_tree_info = f.read().replace(":0.000000;", ";").replace("'", "")
        with open(nwk_tree_outfile, "w") as f:
            f.write(replace_tree_info)

    @staticmethod
    def multifurcate_zero_length_nodes(
        gene_trees_infile: Union[str, Path],
        multifurcate_outfile: Union[str, Path],
        max_tree_num: int = 100,
    ) -> None:
        """Multifurcate zero branch length nodes (For IQ-TREE randomly bifurcated tree)

        Args:
            gene_trees_infile (Union[str, Path]): IQ-TREE gene trees file
            multifurcate_outfile (Union[str, Path]): Multifurcate gene trees file
            max_tree_num (int, optional): Number of max tree output

        Note:
            IQ-TREE randomly bifurcate identical sequences with zero branch length.
            Random bifurcated topology make the DTL reconciliation result worse.
            In order to resolve randomly bifurcated tree using 'treerecs' software,
            multifurcate zero branch length IQ-TREE gene trees to use as treerecs input.
        """
        # Load gene trees
        gene_tree_list: List[EteTree] = []
        with open(gene_trees_infile) as f:
            tree_text_lines = f.read().splitlines()
            for tree_text in tree_text_lines[0:max_tree_num]:
                gene_tree_list.append(EteTree(tree_text))

        # Unroot all zero branch length nodes in gene trees
        gene_tree_text_list = []
        for gene_tree in gene_tree_list:
            gene_tree.set_outgroup(gene_tree.get_midpoint_outgroup())
            for node in gene_tree.traverse(strategy="postorder"):
                descendants_node_list = node.get_descendants()
                if len(descendants_node_list) <= 2:
                    continue
                descendants_total_dist = sum([n.dist for n in descendants_node_list])
                if descendants_total_dist == 0:
                    node.unroot()
            # gene_tree.unroot()
            gene_tree_text_list.append(gene_tree.write(dist_formatter="%1.10f"))

        # Output unrooted gene trees
        with open(multifurcate_outfile, "w") as f:
            f.write("\n".join(gene_tree_text_list))

    @staticmethod
    def unroot_tree(
        gene_trees_infile: Union[str, Path],
        unroot_gene_trees_outfile: Union[str, Path],
    ) -> None:
        """Unroot tree

        Args:
            gene_trees_infile (Union[str, Path]): IQ-TREE gene trees file
            unroot_gene_trees_outfile (Union[str, Path]): Output unroot gene trees file
        """
        # Load gene trees
        gene_tree_list: List[EteTree] = []
        with open(gene_trees_infile) as f:
            tree_text_lines = f.read().splitlines()
            for tree_text in tree_text_lines:
                gene_tree_list.append(EteTree(tree_text))

        # Unroot gene trees
        gene_tree_text_list = []
        for gene_tree in gene_tree_list:
            gene_tree.unroot()
            gene_tree_text_list.append(gene_tree.write(dist_formatter="%1.10f"))

        # Output unrooted gene trees
        with open(unroot_gene_trees_outfile, "w") as f:
            f.write("\n".join(gene_tree_text_list))

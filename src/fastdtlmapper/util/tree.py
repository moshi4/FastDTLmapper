from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Union

from Bio import Phylo, SeqIO
from Bio.Phylo.BaseTree import Tree
from fastdtlmapper.angst.model import NodeEvent


@dataclass
class UtilTree:
    """Tree Utility Class"""

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

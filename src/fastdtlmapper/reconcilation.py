from dataclasses import dataclass
from typing import List

from Bio import Phylo
from Bio.Phylo.BaseTree import Tree

from fastdtlmapper.angst.model import NodeEvent, Transfer
from fastdtlmapper.util import UtilFasta


@dataclass
class Reconciliation:
    """DTL Reconciliation for two or one species Class"""

    @staticmethod
    def two_species(
        species_tree_file: str,
        fasta_file: str,
        los_cost: int = 1,
        trn_cost: int = 3,
    ) -> List[NodeEvent]:
        """DTL reconciliation with two species gene

        Args:
            species_tree_file (str): Species tree file path
            fasta_file (str): Two species gene fasta file path

        Returns:
            List[NodeEvent]: DTL reconciliation node event list
        """
        node_event_list = []
        gene_species_names = UtilFasta(fasta_file).uniq_species_name_list
        tree: Tree = Phylo.read(species_tree_file, "newick")

        # Case: brn & dup in one leaf
        if len(gene_species_names) == 1:
            for node in tree.find_clades():
                if node.name == gene_species_names[0]:
                    node_event_list.append(
                        NodeEvent(node.name, brn_num=1, dup_num=1, gene_num=2)
                    )
                else:
                    node_event_list.append(NodeEvent(node.name, gene_num=0))
            return node_event_list

        # mrca brn & los scenario
        # Most Parsimonious DTL cost calculation
        # TODO: Fix complex code logic more simply
        mrca_node: Tree = tree.common_ancestor(gene_species_names)
        sp1, sp2 = gene_species_names
        los_node_names, no_los_node_names, no_los_no_gene_node_names = [], [], []
        node: Tree
        for node in tree.find_clades():
            if node.name == mrca_node.name:
                continue
            elif not mrca_node.is_parent_of(node.name):
                # No los & no gene node (los_num=0,gene_num=0)
                no_los_no_gene_node_names.append(node.name)
            elif (
                node.name in gene_species_names
                or node.is_parent_of(sp1)
                or node.is_parent_of(sp2)
                or node.is_parent_of(mrca_node.name)
            ):
                # No los node (los_num=0,gene_num=1)
                no_los_node_names.append(node.name)
            elif (
                not node.is_parent_of(sp1)
                and not node.is_parent_of(sp2)
                and node.name not in no_los_no_gene_node_names
            ):
                # Los node (los_num=1,gene_num=0)
                los_node_names.append(node.name)
                no_los_no_gene_node_names.extend(
                    [n.name for n in node.find_clades() if n.name != node.name]
                )

        mrca_brn_los_scenario_cost = los_cost * len(los_node_names)
        if trn_cost >= mrca_brn_los_scenario_cost:
            # Case: mrca brn and other species los gene
            for node in tree.find_clades():
                if node.name == mrca_node.name:
                    node_event_list.append(NodeEvent(node.name, brn_num=1, gene_num=1))
                elif node.name in los_node_names:
                    node_event_list.append(NodeEvent(node.name, los_num=1, gene_num=0))
                elif node.name in no_los_node_names:
                    node_event_list.append(NodeEvent(node.name, los_num=0, gene_num=1))
                elif node.name in no_los_no_gene_node_names:
                    node_event_list.append(NodeEvent(node.name, los_num=0, gene_num=0))
                else:
                    raise ValueError(
                        "ERROR: Two species reconciliation logic may incorrect."
                    )
            return node_event_list
        else:
            # Case: brn in one species and trn other species
            for node in tree.find_clades():
                if node.name == sp1:
                    node_event_list.append(NodeEvent(node.name, brn_num=1, gene_num=1))
                elif node.name == sp2:
                    trn_detail = Transfer(sp1, sp2, direction=False)
                    node_event_list.append(
                        NodeEvent(
                            node.name,
                            trn_num=1,
                            trn_detail_list=[trn_detail],
                            gene_num=1,
                        )
                    )
                else:
                    node_event_list.append(NodeEvent(node.name, gene_num=0))
            return node_event_list

    @staticmethod
    def one_species(
        species_tree_file: str,
        fasta_file: str,
    ) -> List[NodeEvent]:
        """DTL reconciliation with one species gene

        Args:
            species_tree_file (str): Species tree file path
            fasta_file (str): Two species gene fasta file path

        Returns:
            List[NodeEvent]: DTL reconciliation node event list
        """
        node_event_list = []
        gene_species_name = UtilFasta(fasta_file).uniq_species_name_list[0]
        tree: Tree = Phylo.read(species_tree_file, "newick")

        for node in tree.find_clades():
            if node.name == gene_species_name:
                node_event_list.append(NodeEvent(node.name, brn_num=1, gene_num=1))
            else:
                node_event_list.append(NodeEvent(node.name, gene_num=0))

        return node_event_list

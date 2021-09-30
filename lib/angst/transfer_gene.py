from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

from Bio import Phylo
from Bio.Phylo.BaseTree import Tree


@dataclass
class AngstTransferGene:
    """AnGST transfer gene Class"""

    species_tree_file: str
    angst_result_dir: str

    def __post_init__(self):
        self.species_tree_file = Path(self.species_tree_file)
        self.angst_result_dir = Path(self.angst_result_dir)
        self.angst_leaf_file = self.angst_result_dir / "AnGST.leaf"

        self.trn_fromto2gene_id_list = self._parse_leaf()

    def _parse_leaf(self) -> Dict[str, List[str]]:
        """Parse AnGST.leaf file"""
        trn_fromto2gene_id_list = defaultdict(list)

        with open(self.angst_leaf_file) as f:
            line_list = f.read().splitlines()
        tree: Tree = Phylo.read(self.species_tree_file, "newick")
        donor_node, recipient_node = None, None
        for line in line_list:
            if line.startswith("leaf"):
                gene_id = line[5:-1]
            elif line.startswith("	[hgt]:"):
                donor_node_info, recipient_node_info = line[8:].split(" --> ")
                donor_node = tree.common_ancestor(donor_node_info.split("-"))
                recipient_node = tree.common_ancestor(recipient_node_info.split("-"))
            elif line == "":
                if donor_node is not None and recipient_node is not None:
                    trn_fromto = f"{donor_node.name} -> {recipient_node.name}"
                    trn_fromto2gene_id_list[trn_fromto].append(gene_id)
                donor_node, recipient_node = None, None

        return trn_fromto2gene_id_list

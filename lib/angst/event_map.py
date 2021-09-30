from dataclasses import dataclass
from pathlib import Path
from typing import Dict

from Bio import Phylo
from Bio.Phylo.BaseTree import Tree
from lib.angst.model import NodeEvent, Transfer


@dataclass
class AngstEventMap:
    """AnGST Event Map Class"""

    species_tree_file: str
    angst_result_dir: str

    def __post_init__(self):
        self.species_tree_file = Path(self.species_tree_file)
        self.angst_result_dir = Path(self.angst_result_dir)
        self.angst_event_file = self.angst_result_dir / "AnGST.events"

        self.nodeid2node_event = self._parse_events()

    def _parse_events(self) -> Dict[str, NodeEvent]:
        """Parse AnGST.events file"""
        # Get initialized node event dict
        nodeid2node_event = self._get_blank_node_event_dict()

        # Parse AnGST event file & Map event
        with open(self.angst_event_file) as f:
            event_line_list = f.read().splitlines()
        tree: Tree = Phylo.read(self.species_tree_file, "newick")
        for event_line in event_line_list:
            event_type, event_node_info = event_line[1:4], event_line[7:]

            if event_type == "hgt":
                donor_node_info, recipient_node_info = event_node_info.split(" --> ")
                donor_node = tree.common_ancestor(donor_node_info.split("-"))
                recipient_node = tree.common_ancestor(recipient_node_info.split("-"))
                nodeid2node_event[recipient_node.name].trn_num += 1
                nodeid2node_event[recipient_node.name].trn_detail_list.append(
                    Transfer(donor_node.name, recipient_node.name)
                )
                continue

            lca_node = tree.common_ancestor(event_node_info.split("-"))
            if event_type == "brn":
                nodeid2node_event[lca_node.name].brn_num += 1
            elif event_type == "dup":
                nodeid2node_event[lca_node.name].dup_num += 1
            elif event_type == "los":
                nodeid2node_event[lca_node.name].los_num += 1

        # Countup ancestor gene num from event map info
        for node_event in nodeid2node_event.values():
            root_to_countup_node_path = [
                tree.root,
                *tree.get_path(node_event.node_id),
            ]
            gene_count = 0
            for node in root_to_countup_node_path:
                node_event = nodeid2node_event[node.name]
                gene_count += node_event.gain_num - node_event.los_num
            node_event.gene_num = gene_count

        return nodeid2node_event

    def _get_blank_node_event_dict(self) -> Dict[str, NodeEvent]:
        """Get blank node event dict for event map initialization"""
        tree: Tree = Phylo.read(self.species_tree_file, "newick")
        nodeid2node_event = {}
        for node in tree.find_clades():
            nodeid2node_event[node.name] = NodeEvent(node_id=node.name)
        return nodeid2node_event

    def write_tree(self, outfile: str, map_type: str) -> None:
        """Write event map tree newick file

        Args:
            outfile (str): Output tree file path
            map_type (str): mapping type ("gain-loss" or "dtl")
        """
        tree: Tree = Phylo.read(self.species_tree_file, "newick")
        for node in tree.find_clades():
            node_event = self.nodeid2node_event[node.name]
            if map_type == "gain-loss":
                node.name = node.name + " | " + node_event.as_gain_loss_text
            elif map_type == "dtl":
                node.name = node.name + " | " + node_event.as_dtl_text
            else:
                raise ValueError("type must be 'gain-loss' or 'dtl'")

        Phylo.write(tree, outfile, "newick")
        with open(outfile) as f:
            replace_tree_info = f.read().replace(":0.00000;", ";").replace("'", "")
        with open(outfile, "w") as f:
            f.write(replace_tree_info)

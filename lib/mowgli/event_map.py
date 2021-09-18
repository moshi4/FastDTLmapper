from dataclasses import dataclass
from typing import List

from Bio import Phylo
from lib.mowgli.model import NodeEvent, Transfer
from lib.mowgli.tree import MowgliTree
from lib.mowgli.xml_parser import MowgliXmlParser


@dataclass
class MowgliEventMap:
    """Mowgli Event Map Class"""

    mowgli_species_tree_file: str
    mowgli_xml_file: str

    def __post_init__(self):
        self.mowgli_xml = MowgliXmlParser(self.mowgli_xml_file)
        self.mowgli_tree = MowgliTree(self.mowgli_species_tree_file)

    def write_tree(self, outfile: str, map_type: str) -> None:
        """Write event map tree newick file

        Args:
            outfile (str): Output tree file path
            map_type (str): mapping type ("gain-loss" or "dtl")
        """
        named_tree = self.mowgli_tree.named_tree
        for node in named_tree.find_clades():
            node_event = self.get_node_event(node.name)
            if map_type == "gain-loss":
                node.name = node.name + " | " + node_event.as_gain_loss_text
            elif map_type == "dtl":
                node.name = node.name + " | " + node_event.as_dtl_text
            else:
                raise ValueError("type must be 'gain-loss' or 'dtl'")

        Phylo.write(named_tree, outfile, "newick")
        with open(outfile) as f:
            replace_tree_info = f.read().replace(":0.00000;", ";").replace("'", "")
        with open(outfile, "w") as f:
            f.write(replace_tree_info)

    def get_all_node_event(self) -> List[NodeEvent]:
        """Get all node event

        Returns:
            List[NodeEvent]: All node event
        """
        node_event_list = []
        for node_name in self.mowgli_tree.node_num2name.values():
            node_event_list.append(self.get_node_event(node_name))
        return node_event_list

    def get_node_event(self, target_node_id: str) -> NodeEvent:
        """Get target node event

        Args:
            target_node_id (str): Target node id

        Returns:
            NodeEvent: Target node event
        """
        # Get DTL node event
        node_name2num = {
            name: int(num) for num, name in self.mowgli_tree.node_num2name.items()
        }
        node_event = self.mowgli_xml.get_node_event(node_name2num[target_node_id])

        # Convert transfer detail 'node num' to 'node name'
        trn_detail_list = []
        for trn in node_event.trn_detail_list:
            node_num2name = self.mowgli_tree.node_num2name
            donor_node_id = node_num2name[str(trn.donor_node_id)]
            recipient_node_id = node_num2name[str(trn.recipient_node_id)]
            trn_detail_list.append(Transfer(donor_node_id, recipient_node_id))

        # Count node gene number
        named_tree = self.mowgli_tree.named_tree
        root_to_target_node_path = [
            named_tree.root,
            *named_tree.get_path(target_node_id),
        ]
        gene_count = 0
        for node in root_to_target_node_path:
            node_event = self.mowgli_xml.get_node_event(node_name2num[node.name])
            gene_count += node_event.gain_num - node_event.los_num

        # Update node event contents
        node_event.node_id = target_node_id
        node_event.gene_num = gene_count
        node_event.trn_detail_list = trn_detail_list

        return node_event

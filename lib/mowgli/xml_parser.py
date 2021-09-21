from dataclasses import dataclass
from pathlib import Path
from typing import List, Union

from bs4 import BeautifulSoup
from lib.mowgli.model import NodeEvent, Transfer


@dataclass
class MowgliXmlParser:
    """Mowgli Xml Class"""

    mowgli_xml_file: Union[str, Path]

    def __post_init__(self):
        self._soup = self._read_xml()
        self._set_params()
        self._brn_node_id = self._get_brn_node_id()
        self._dup_node_id_list = self._get_dup_node_id_list()
        self._los_node_id_list = self._get_los_node_id_list()
        self._trn_list = self._get_trn_list()

    def get_node_event(self, target_node_id: int) -> NodeEvent:
        """Get target node DTL event

        Args:
            target_node_id (int): Target node id

        Returns:
            NodeEvent: Target node event
        """
        dup_num = self._dup_node_id_list.count(target_node_id)
        los_num = self._los_node_id_list.count(target_node_id)
        trn_num = 0
        trn_detail_list = []
        for trn in self._trn_list:
            if target_node_id == trn.recipient_node_id:
                trn_num += 1
                trn_detail_list.append(trn)
        brn_num = 1 if self._brn_node_id == target_node_id else 0
        return NodeEvent(
            target_node_id, brn_num, dup_num, los_num, trn_num, trn_detail_list
        )

    def _read_xml(self) -> BeautifulSoup:
        """Read Mowgli xml result

        Returns:
            BeautifulSoup: Beautifulsoup xml parsed object
        """
        with open(self.mowgli_xml_file) as f:
            # Replace invalid parsable character ':' to '_'
            xml_string = f.read().replace(":", "_")
        return BeautifulSoup(xml_string, "lxml")

    def _set_params(self) -> None:
        """Set cost paramergers"""
        param_value_list = [param["value"] for param in self._soup.select("rec_param")]
        self.dup_cost = int(param_value_list[0])
        self.los_cost = int(param_value_list[1])
        self.spc_cost = int(param_value_list[2])
        self.trn_cost = int(param_value_list[3])
        self.total_cost = int(self._soup.select_one("rec_totalCost").text)

    def _get_brn_node_id(self) -> int:
        """Get gene born node id

        Returns:
            int: Gene born node id
        """
        normal_brn_node_tag = self._soup.select_one("rec_locationSp")
        if normal_brn_node_tag is not None:
            return int(normal_brn_node_tag.text)
        else:
            # Case: Transfer event only
            transfer_only_brn_node_tag = self._soup.select_one("rec_originSp")
            return int(transfer_only_brn_node_tag.text)

    def _get_dup_node_id_list(self) -> List[int]:
        """Get duplication node id list

        Returns:
            List[int]: Duplication node id list
        """
        dup_node_id_list = []
        for dup_record in self._soup.select("rec_duplication"):
            dup_node_id = int(dup_record.select_one("rec_locationSp").text)
            dup_node_id_list.append(dup_node_id)
        return dup_node_id_list

    def _get_los_node_id_list(self) -> List[int]:
        """Get loss node id list

        Returns:
            List[int]: Loss node id list
        """
        los_node_id_list = []
        for los_record in self._soup.select("rec_loss"):
            los_node_id = int(los_record.select_one("rec_locationSp").text)
            los_node_id_list.append(los_node_id)
        return los_node_id_list

    def _get_trn_list(self) -> List[Transfer]:
        """Get transfer donor & recipient node id list

        Returns:
            List[Transfer]: Transfer node id list
        """
        trn_node_id_list = []
        for trn_record in self._soup.select("rec_transfer"):
            trn_donor_node_id = int(trn_record.select_one("rec_originSp").text)
            trn_recipient_node_id = int(trn_record.select_one("rec_recipientSp").text)
            trn_node_id_list.append(Transfer(trn_donor_node_id, trn_recipient_node_id))
        return trn_node_id_list

    def __str__(self) -> str:
        ret_val = ""
        ret_val += f"dup_cost={self.dup_cost}\n"
        ret_val += f"los_cost={self.los_cost}\n"
        ret_val += f"spc_cost={self.spc_cost}\n"
        ret_val += f"trn_cost={self.trn_cost}\n"
        ret_val += f"total_cost={self.total_cost}\n"
        ret_val += f"brn_node_id={self._brn_node_id}\n"
        ret_val += f"dup_node_id_list={self._dup_node_id_list}\n"
        ret_val += f"los_node_id_list={self._los_node_id_list}\n"
        ret_val += f"trn_list={self._trn_list}\n"
        return ret_val

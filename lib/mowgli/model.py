from dataclasses import dataclass
from typing import List, Optional, Union


@dataclass
class Transfer:
    """Transfer DataClass"""

    donor_node_id: Union[int, str]
    recipient_node_id: Union[int, str]

    @property
    def as_text(self) -> str:
        """Transfer as text format (e.g. 1 -> 2)"""
        return f"{self.donor_node_id} -> {self.recipient_node_id}"


@dataclass
class NodeEvent:
    """Node Event DataClass"""

    node_id: int
    dup_num: int
    los_num: int
    trn_num: int
    trn_detail_list: List[Transfer]
    brn_num: bool
    gene_num: Optional[int] = None

    @property
    def gain_num(self) -> int:
        """Number of (Brn + Dup + Trn) event"""
        return self.brn_num + self.dup_num + self.trn_num

    @property
    def as_dtl_text(self) -> str:
        """DTL event text"""
        text = ""
        text += f"{self.gene_num} ["
        text += f"brn={self.brn_num} "
        text += f"dup={self.dup_num} "
        text += f"los={self.los_num} "
        text += f"trn={self.trn_num}]"
        return text

    @property
    def as_gain_loss_text(self) -> str:
        """Gain-Loss event text"""
        text = ""
        text += f"{self.gene_num} ["
        text += f"gain={self.gain_num} "
        text += f"los={self.los_num}]"
        return text
